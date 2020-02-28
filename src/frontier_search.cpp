#include <explore/frontier_search.h>

#include <mutex>
#include <math.h>
#include <ros/ros.h>
#include <costmap_2d/cost_values.h>
#include <costmap_2d/costmap_2d.h>
#include <geometry_msgs/Point.h>
#include <tf/transform_datatypes.h>
#include <explore/costmap_tools.h>


namespace frontier_exploration
{
using costmap_2d::LETHAL_OBSTACLE;
using costmap_2d::NO_INFORMATION;
using costmap_2d::FREE_SPACE;

FrontierSearch::FrontierSearch(costmap_2d::Costmap2D* costmap,
                               double potential_scale, double gain_scale,
                               double min_frontier_size,
                               double min_x, double min_y, double max_x,
                               double max_y, geometry_msgs::PoseStamped initial_pose)
  : costmap_(costmap)
  , potential_scale_(potential_scale)
  , gain_scale_(gain_scale)
  , min_frontier_size_(min_frontier_size)
  , min_x_(min_x)
  , min_y_(min_y)
  , max_x_(max_x)
  , max_y_(max_y)
  , initial_pose_(initial_pose)
  , roll_(0.0)
  , pitch_(0.0)
  , theta_(0.0)
{
  tf::Pose pose;
  tf::poseMsgToTF(initial_pose_.pose, pose);
  theta_ = tf::getYaw(pose.getRotation());

}

std::vector<Frontier> FrontierSearch::searchFrom(geometry_msgs::Point position)
{
  std::vector<Frontier> frontier_list;

  // Sanity check that robot is inside costmap bounds before searching
  unsigned int mx, my;
  if (!costmap_->worldToMap(position.x, position.y, mx, my)) {
    ROS_ERROR("Robot out of costmap bounds, cannot search for frontiers");
    return frontier_list;
  }

  // make sure map is consistent and locked for duration of search
  std::lock_guard<costmap_2d::Costmap2D::mutex_t> lock(*(costmap_->getMutex()));

  map_ = costmap_->getCharMap();
  size_x_ = costmap_->getSizeInCellsX();
  size_y_ = costmap_->getSizeInCellsY();

  // initialize flag arrays to keep track of visited and frontier cells
  std::vector<bool> frontier_flag(size_x_ * size_y_, false);
  std::vector<bool> visited_flag(size_x_ * size_y_, false);

  // initialize breadth first search
  std::queue<unsigned int> bfs;

  // find closest clear cell to start search
  unsigned int clear, pos = costmap_->getIndex(mx, my);
  if (nearestCell(clear, pos, FREE_SPACE, *costmap_)) {
    bfs.push(clear);
  } else {
    bfs.push(pos);
    ROS_WARN("Could not find nearby clear cell to start search");
  }
  visited_flag[bfs.front()] = true;

  while (!bfs.empty()) {
    unsigned int idx = bfs.front();
    bfs.pop();

    // iterate over 4-connected neighbourhood
    for (unsigned nbr : nhood4(idx, *costmap_)) {
      // add to queue all free, unvisited cells, use descending search in case
      // initialized on non-free cell
      if (map_[nbr] <= map_[idx] && !visited_flag[nbr]) {
        visited_flag[nbr] = true;
        bfs.push(nbr);
        // check if cell is new frontier cell (unvisited, NO_INFORMATION, free
        // neighbour)
      } else if (isNewFrontierCell(nbr, frontier_flag)) {
        frontier_flag[nbr] = true;
        Frontier new_frontier = buildNewFrontier(nbr, pos, frontier_flag);
        if (new_frontier.size * costmap_->getResolution() >=
            min_frontier_size_) {
          frontier_list.push_back(new_frontier);
        }
      }
    }
  }

  Boundary boundary = getMapBoundary(min_x_, min_y_, max_x_, max_y_, 
                                     initial_pose_.pose.position.x, 
                                     initial_pose_.pose.position.y, 
                                     theta_);

  //now that we have our basic frontier list, we can prune the frontiers that
  // are not within the desired boundaries
  getBoundedFrontierList( frontier_list, boundary);

  // set costs of frontiers
  for (auto& frontier : frontier_list) {
    frontier.cost = frontierCost(frontier);
  }
  std::sort(
      frontier_list.begin(), frontier_list.end(),
      [](const Frontier& f1, const Frontier& f2) { return f1.cost < f2.cost; });

  return frontier_list;
}

Frontier FrontierSearch::buildNewFrontier(unsigned int initial_cell,
                                          unsigned int reference,
                                          std::vector<bool>& frontier_flag)
{
  // initialize frontier structure
  Frontier output;
  output.centroid.x = 0;
  output.centroid.y = 0;
  output.size = 1;
  output.min_distance = std::numeric_limits<double>::infinity();

  // record initial contact point for frontier
  unsigned int ix, iy;
  costmap_->indexToCells(initial_cell, ix, iy);
  costmap_->mapToWorld(ix, iy, output.initial.x, output.initial.y);

  // push initial gridcell onto queue
  std::queue<unsigned int> bfs;
  bfs.push(initial_cell);

  // cache reference position in world coords
  unsigned int rx, ry;
  double reference_x, reference_y;
  costmap_->indexToCells(reference, rx, ry);
  costmap_->mapToWorld(rx, ry, reference_x, reference_y);

  while (!bfs.empty()) {
    unsigned int idx = bfs.front();
    bfs.pop();

    // try adding cells in 8-connected neighborhood to frontier
    for (unsigned int nbr : nhood8(idx, *costmap_)) {
      // check if neighbour is a potential frontier cell
      if (isNewFrontierCell(nbr, frontier_flag)) {
        // mark cell as frontier
        frontier_flag[nbr] = true;
        unsigned int mx, my;
        double wx, wy;
        costmap_->indexToCells(nbr, mx, my);
        costmap_->mapToWorld(mx, my, wx, wy);

        geometry_msgs::Point point;
        point.x = wx;
        point.y = wy;
        output.points.push_back(point);

        // update frontier size
        output.size++;

        // update centroid of frontier
        output.centroid.x += wx;
        output.centroid.y += wy;

        // determine frontier's distance from robot, going by closest gridcell
        // to robot
        double distance = sqrt(pow((double(reference_x) - double(wx)), 2.0) +
                               pow((double(reference_y) - double(wy)), 2.0));
        if (distance < output.min_distance) {
          output.min_distance = distance;
          output.middle.x = wx;
          output.middle.y = wy;
        }

        // add to queue for breadth first search
        bfs.push(nbr);
      }
    }
  }

  // average out frontier centroid
  output.centroid.x /= output.size;
  output.centroid.y /= output.size;
  return output;
}

bool FrontierSearch::isNewFrontierCell(unsigned int idx,
                                       const std::vector<bool>& frontier_flag)
{
  // check that cell is unknown and not already marked as frontier
  if (map_[idx] != NO_INFORMATION || frontier_flag[idx]) {
    return false;
  }

  // frontier cells should have at least one cell in 4-connected neighbourhood
  // that is free
  for (unsigned int nbr : nhood4(idx, *costmap_)) {
    if (map_[nbr] == FREE_SPACE) {
      return true;
    }
  }

  return false;
}

Boundary FrontierSearch::getMapBoundary(const double min_x, const double min_y, 
                                          const double max_x, const double max_y, 
                                          const double init_x, const double init_y, 
                                          const double theta)
{

  Boundary boundary;

  boundary.min_x = init_x + cos(theta) * min_x + sin(theta) * min_y;
  boundary.min_y = init_y + sin(theta) * min_x - cos(theta) * min_y;

  boundary.max_x = init_x + cos(theta) * max_x + sin(theta) * max_y;
  boundary.max_y = init_y + sin(theta) * max_x - cos(theta) * max_y;

  if (boundary.max_x <= boundary.min_x)
  {
    double temp_x = boundary.max_x;
    boundary.max_x = boundary.min_x;
    boundary.min_x = temp_x;
  }

  if (boundary.max_y <= boundary.min_y)
  {
    double temp_y = boundary.max_y;
    boundary.max_y = boundary.min_y;
    boundary.min_y = temp_y;
  }


  return boundary;

}

void FrontierSearch::getBoundedFrontierList(std::vector<Frontier> &frontier_list,
                                                         const Boundary boundary)
{

  std::cout << "Boundaries: Xmin: " << boundary.min_x << " Xmax:" << boundary.max_x << " Ymin: " << boundary.min_y << " Ymax:" << boundary.max_y << " \n";

  for (int i = 0; i < frontier_list.size(); i++)  
  {

    std::cout << "Frontier at X: " << frontier_list[i].centroid.x << " Y:" << frontier_list[i].centroid.y << " \n";

    if (frontier_list[i].centroid.x <= boundary.min_x || 
        frontier_list[i].centroid.x >= boundary.max_x || 
        frontier_list[i].centroid.y <= boundary.min_y ||
        frontier_list[i].centroid.y >= boundary.max_y )
    {
      std::cout << "...deleted \n";
      frontier_list.erase(frontier_list.begin() + i);
    }
    else
    {
      std::cout << " ...not deleted\n";
    }
  }

}

double FrontierSearch::frontierCost(const Frontier& frontier)
{
  return (potential_scale_ * (1.0/frontier.min_distance) *
          costmap_->getResolution()) -
         (gain_scale_ * frontier.size * costmap_->getResolution());
}

}
