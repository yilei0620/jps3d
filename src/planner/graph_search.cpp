#include <jps3d/planner/graph_search.h>

using namespace JPS;

GraphSearch::GraphSearch(const char* cMap, int xDim, int yDim,  double eps, bool verbose) :
  cMap_(cMap), xDim_(xDim), eps_(eps), verbose_(verbose) 
{ 
  hm_.reserve(xDim_ * yDim_);
  seen_.resize(xDim_ * yDim_, false);
}


inline int GraphSearch::coordToId(int x, int y) const {
  return x + y*xDim_;
}

inline bool GraphSearch::isFree(int x, int y) const {
  return x >= 0 && x < xDim_ && y >= 0 && y < yDim_ &&
    cMap_[coordToId(x, y)] == val_free_;
}

inline double GraphSearch::get_heur(int x, int y) const {
  return eps_ * std::sqrt((x - xGoal_) * (x - xGoal_) + (y - yGoal_) * (y - yGoal_));
}

bool GraphSearch::plan(int xStart, int yStart, int xGoal, int yGoal, int max_expand) 
{
  pq_.clear();
  path_.clear();
  hm_.reserve(xDim_ * yDim_);
  seen_.resize(xDim_ * yDim_, false);

  // Set goal
  int start_id = coordToId(xStart, yStart);
  int goal_id = coordToId(xGoal, yGoal);
  xGoal_ = xGoal; yGoal_ = yGoal;

  // Initialize start node
  StatePtr currNode_ptr = std::make_shared<State>(State(coordToId(xStart, yStart), xStart, yStart, 0, 0));
  currNode_ptr->g = 0;
  currNode_ptr->h = get_heur(xStart, yStart);
  currNode_ptr->heapkey = pq_.push(currNode_ptr);
  currNode_ptr->opened = true;
  hm_[currNode_ptr->id] = currNode_ptr;
  seen_[currNode_ptr->id] = true;
 
  int expand_iteration = 0;
  while(true)
  {
    expand_iteration++;
    // get element with smallest cost
    currNode_ptr = pq_.top(); pq_.pop();
    currNode_ptr->closed = true; // Add to closed list

    if(currNode_ptr->id == goal_id) {
      if(verbose_)
        printf("Goal Reached!!!!!!\n\n");
      break;
    }

    std::vector<int> succ_ids;
    std::vector<double> succ_costs;
    // Get successors
    getSucc(currNode_ptr, succ_ids, succ_costs);

    // Process successors
    for( int s = 0; s < (int) succ_ids.size(); s++ )
    {
      //see if we can improve the value of succstate
      StatePtr& child_ptr = hm_[succ_ids[s]];
      double tentative_gval = currNode_ptr->g + succ_costs[s];

      if( tentative_gval < child_ptr->g )
      {
        child_ptr->parentId = currNode_ptr->id;  // Assign new parent
        child_ptr->g = tentative_gval;    // Update gval

        //double fval = child_ptr->g + child_ptr->h;

        // if currently in OPEN, update
        if( child_ptr->opened && !child_ptr->closed)
          pq_.increase( child_ptr->heapkey );       // update heap
        // if currently in CLOSED
        else if( child_ptr->opened && child_ptr->closed)
        {
          printf("ASTAR ERROR!\n");
        }
        else // new node, add to heap
        {
          child_ptr->heapkey = pq_.push(child_ptr);
          child_ptr->opened = true;
        }
      } //
    } // Process successors    


    if(max_expand > 0 && expand_iteration >= max_expand) {
      printf("MaxExpandStep [%d] Reached!!!!!!\n\n", max_expand);
      return false;
    }

    if( pq_.empty()) {
      printf("Priority queue is empty!!!!!!\n\n");
      return false;
    }
  }

  if(verbose_) {
    printf("goal g: %f, h: %f!\n", currNode_ptr->g, currNode_ptr->h);
    printf("Expand [%d] nodes!\n", expand_iteration);
  }

  path_ = recoverPath(currNode_ptr, start_id);

  return true;
}

std::vector<StatePtr> GraphSearch::recoverPath(StatePtr node, int start_id) {
  std::vector<StatePtr> path;
  path.push_back(node);
  while(node && node->id != start_id) {
    node = hm_[node->parentId];
    path.push_back(node);
  }

  return path;
}

void GraphSearch::getSucc(const StatePtr& curr, std::vector<int>& succ_ids, std::vector<double>& succ_costs) {
  for(int dx = -1; dx <= 1; dx ++) {
    for(int dy = -1; dy <= 1; dy ++) {
      if(dx == 0 && dy == 0)
        continue;
      int new_x = curr->x + dx;
      int new_y = curr->y + dy;
      if(!isFree(new_x, new_y))
        continue;

      int new_id = coordToId(new_x, new_y);
      if(!seen_[new_id]) {
        seen_[new_id] = true;
        hm_[new_id] = std::make_shared<State>(new_id, new_x, new_y, dx, dy);
        hm_[new_id]->h = get_heur(new_x, new_y);
      }

      succ_ids.push_back(new_id);
      succ_costs.push_back(std::sqrt(dx*dx+dy*dy));
    }
  }
}

std::vector<StatePtr> GraphSearch::getPath() const {
  return path_;
}
