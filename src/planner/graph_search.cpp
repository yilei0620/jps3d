#include <jps3d/planner/graph_search.h>

using namespace JPS;

GraphSearch::GraphSearch(const char* cMap, int xDim, int yDim,  double eps, bool verbose) :
  cMap_(cMap), xDim_(xDim), yDim_(yDim), eps_(eps), verbose_(verbose) 
{ 
  hm_.resize(xDim_ * yDim_);
  seen_.resize(xDim_ * yDim_, false);
}


inline int GraphSearch::coordToId(int x, int y) const {
  return x + y*xDim_;
}

inline bool GraphSearch::isFree(int x, int y) const {
  return x >= 0 && x < xDim_ && y >= 0 && y < yDim_ &&
    cMap_[coordToId(x, y)] == val_free_;
}

inline bool GraphSearch::isOccupied(int x, int y) const {
  return x >= 0 && x < xDim_ && y >= 0 && y < yDim_ &&
    cMap_[coordToId(x, y)] > val_free_;
}



inline double GraphSearch::get_heur(int x, int y) const {
  return eps_ * std::sqrt((x - xGoal_) * (x - xGoal_) + (y - yGoal_) * (y - yGoal_));
}

bool GraphSearch::plan(int xStart, int yStart, int xGoal, int yGoal, int max_expand) 
{
  pq_.clear();
  path_.clear();
  hm_.resize(xDim_ * yDim_);
  seen_.resize(xDim_ * yDim_, false);

  // Set goal
  int start_id = coordToId(xStart, yStart);
  int goal_id = coordToId(xGoal, yGoal);
  xGoal_ = xGoal; yGoal_ = yGoal;

  // Initialize start node
  StatePtr currNode_ptr = std::make_shared<State>(State(start_id, xStart, yStart, 0, 0));
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

    //printf("expand: %d, %d\n", currNode_ptr->x, currNode_ptr->y);
    std::vector<int> succ_ids;
    std::vector<double> succ_costs;
    // Get successors
    if(!jps_)
      getSucc(currNode_ptr, succ_ids, succ_costs);
    else
      getJpsSucc(currNode_ptr, succ_ids, succ_costs);

    //printf("size of succs: %zu\n", succ_ids.size());
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
          //printf("add to open set: %d, %d\n", child_ptr->x, child_ptr->y);
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

void GraphSearch::getJpsSucc(const StatePtr& curr, std::vector<int>& succ_ids, std::vector<double>& succ_costs) {
  std::vector<std::array<int, 2>> ns = prune(curr->x, curr->y, curr->dx, curr->dy);

  for(const auto& it: ns) {
    int new_x, new_y;
    if(jump(curr->x, curr->y, it[0], it[1], new_x, new_y)) {
      int new_id = coordToId(new_x, new_y);
      if(!seen_[new_id]) {
        seen_[new_id] = true;
        hm_[new_id] = std::make_shared<State>(new_id, new_x, new_y, it[0], it[1]); 
        hm_[new_id]->h = get_heur(new_x, new_y);
      }

      succ_ids.push_back(new_id);
      succ_costs.push_back(std::sqrt((new_x - curr->x) * (new_x - curr->x) + 
            (new_y - curr->y) * (new_y  - curr->y)));
    }
  }
}


bool GraphSearch::jump(int x, int y, int dx, int dy, int& new_x, int& new_y ) {
  new_x = x + dx;
  new_y = y + dy;
  if (!isFree(new_x, new_y))
    return false;

  //printf("dx: %d, dy: %d\n", dx, dy);
  if (new_x ==  xGoal_ && new_y == yGoal_)
    return true;

  if (hasForced(new_x, new_y, dx, dy))
    return true;

  // Diagonal jump
  if (std::abs(dx) + std::abs(dy) == 2) {
    int i = 3 * (dx + 1) + (dy + 1);
    for (const auto &it : add_map_[i]) {
      int new_new_x, new_new_y;
      if (!(it[0] == dx && it[1] == dy) && jump(new_x, new_y, it[0], it[1], new_new_x, new_new_y)) 
        return true;
    }
  }

  return jump(new_x, new_y, dx, dy, new_x, new_y);
}

bool GraphSearch::hasForced(int x, int y, int dx, int dy) {
  int i = 3*(dx + 1) + (dy + 1);
  for (const auto &it : obs_map_[i]) {
    if (isOccupied(x + it[0], y + it[1])) {
      return true;
    }
  }

  return false;
}


std::vector<std::array<int, 2>> GraphSearch::prune(int x, int y, int dx, int dy) {
  std::vector<std::array<int, 2>> ns;
  if (dx == 0 && dy == 0) {
    for(int i = -1; i <= 1; i++) {
      for(int j = -1; j <= 1; j++) {
	if(i == 0 && j == 0)
	  continue;
	ns.push_back(std::array<int, 2>{i, j});
      }
    }
    return ns;
  }
  int dir = 3 * (dx + 1) + (dy + 1);
  for (const auto &it : add_map_[dir]) 
    ns.push_back(it);

  if (std::abs(dx) + std::abs(dy) == 1) {
    for (const auto &it : obs_map_[dir]) {
      if (isOccupied(x + it[0], y + it[1]) && isFree(x + it[0] + dx, y + it[1] + dy)) 
        ns.push_back(std::array<int ,2> {it[0] + dx, it[1] + dy});
    }
  }

  else if (std::abs(dx) + std::abs(dy) == 2) {
    for (const auto &it : obs_map_[dir]) {
      if (isOccupied(x + it[0], y + it[1]) && isFree(x + 2 * it[0] + dx, y + 2 * it[1] + dy))
        ns.push_back(std::array<int, 2> {2 * it[0] + dx, 2 * it[1] + dy});
    }
  }

  return ns;
}

std::vector<StatePtr> GraphSearch::getPath() const {
  return path_;
}

std::vector<StatePtr> GraphSearch::getOpenedState() const {
  std::vector<StatePtr> ss;
  for(const auto& it: hm_) {
    if(it && it->opened)
      ss.push_back(it);
  }
  return ss;
}

void GraphSearch::useJps() {
  add_map_.resize(9);
  obs_map_.resize(9);
  int cnt = 0;
  //**** neightbors and moves
  for(int x = -1; x <= 1; x++) {
    for(int y = -1; y <= 1; y++) {
      std::vector<std::array<int, 2>> add;
      std::vector<std::array<int, 2>> obs;

      //*** straight move
      if ((std::abs(x) == 1 && std::abs(y) != 1) ||
	  (std::abs(x) != 1 && std::abs(y) == 1)) {
	add.push_back(std::array<int, 2>{x, y}); 

	for(int dx = -1; dx <= 1; dx++) {
	  for(int dy = -1; dy <= 1; dy++) {
	    if (x * dx + y * dy == 0 && !(dx == 0 && dy == 0))
	      obs.push_back(std::array<int, 2>{dx, dy});
	  }   
	}
      }

      //*** one-d diagonal move
      else if (std::abs(x) + std::abs(y) == 2) {
	add.push_back(std::array<int, 2>{x, y});
	add.push_back(std::array<int, 2>{x, 0});
	add.push_back(std::array<int, 2>{0, y});
	obs.push_back(std::array<int, 2>{-x, 0});
	obs.push_back(std::array<int, 2>{0, -y});
      }

      //printf("dx: %d, dy: %d --- \n", x, y);
      add_map_[cnt] = add;
      //for(const auto& it: add)
	//printf("add: %d, %d\n", it[0], it[1]);
      obs_map_[cnt] = obs;
      //for(const auto& it: obs)
	//printf("obs: %d, %d\n", it[0], it[1]);

      cnt ++;
    }
  }
  
  jps_ = true;

  if(verbose_)
    printf("Use JPS!\n\n");
}
