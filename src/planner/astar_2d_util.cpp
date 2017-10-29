#include <jps3d/planner/astar_2d_util.h>

using namespace JPS;

AStar2DUtil::AStar2DUtil(bool verbose) {
  _planner_verbose = verbose;
  if(_planner_verbose)
    printf(ANSI_COLOR_CYAN "PLANNER VERBOSE ON\n" ANSI_COLOR_RESET);

  _status = 0;
  // -------
  // 0, normal
  // -1, goal is occupied or cannot find a path
  // 1, too close, dont plan
}

bool AStar2DUtil::plan(const Vec3f &start, const Vec3f &goal, decimal_t eps) {
  if(_planner_verbose){
    std::cout <<"Start: " << start.transpose() << std::endl;
    std::cout <<"Goal:  " << goal.transpose() << std::endl;
    std::cout <<"Epsilon:  " << eps << std::endl;
  }

  _path.clear();
  _raw_path.clear();
  _status = 0;

  _start_int = _map_util->floatToInt(start);

  if (!_map_util->isFree(_start_int)) {
    if(_planner_verbose)
      printf(ANSI_COLOR_YELLOW "start is not free! cleared" ANSI_COLOR_RESET "\n");
    if(_map_util->isOutSide(_start_int)) {
      Vec3i dim = _map_util->getDim();
      if(_planner_verbose)
        printf(ANSI_COLOR_RED "Start is outside [%d, %d, %d], dim: [%d, %d, %d]" ANSI_COLOR_RESET "\n",
               _start_int(0), _start_int(1), _start_int(2),
               dim(0), dim(1), dim(2));
      _status = 1;
      return false;
    }
    //_map_util->clearAround(_start_int);
    _status = 1;
    return false;
  }


  _goal_int = _map_util->floatToInt(goal);
  _goal_outside = false;
  if((_goal_int - _start_int).topRows<2>().norm() == 0)
  {
    if(_planner_verbose)
      printf(ANSI_COLOR_RED "start and goal are too close! " ANSI_COLOR_RESET "\n");
    _status = 1;
    return false;
  }

  if (_map_util->isOutSideXYZ(_goal_int, 0) ||
      _map_util->isOutSideXYZ(_goal_int, 1))
    _goal_outside = true;
  else
  {
    if(_map_util->isOutSideXYZ(_goal_int, 2)) {
      _goal_int(2) = _start_int(2);
      if(_planner_verbose)
        printf(ANSI_COLOR_CYAN "change goal z! " ANSI_COLOR_RESET "\n");
    }

    if (!_map_util->isFree(_goal_int) ) {
      Vec3f g = _map_util->intToFloat(_goal_int);
      vec_Vec3i poses = _map_util->rayTrace(start, g);
      int id = -1;
      for(int i = poses.size() - 1; i > 0; i-- ){
        if(_map_util->isFree(poses[i]) &&
            _map_util->isFree(poses[i-1]) )
        {
          id = i - 1;
          break;
        }
      }
      if(id >= 0) {
        _goal_int = poses[id];
        if(_planner_verbose)
          printf(ANSI_COLOR_CYAN "change goal! " ANSI_COLOR_RESET "\n");

      }
      else{
        if(_planner_verbose)
          printf(ANSI_COLOR_RED "goal is not free! Abort! " ANSI_COLOR_RESET "\n");
        _status = -1;
        return false;
      }
    }
  }

  Vec3i dim = _map_util->getDim();

  std::vector<char> cmap(dim(0)*dim(1));
  for( int y = 0; y < dim(1); ++y)
    for( int x = 0; x < dim(0); ++x)
    {
      if(_map_util->isOccupied(Vec3i(x,y,0)))
        cmap[x +y*dim(0)] = 1;
      else
        cmap[x +y*dim(0)] = 0;
    }

  JPS::GraphSearch gs(cmap.data(), dim(0), dim(1), eps, _planner_verbose);

  gs.plan(_start_int(0), _start_int(1), _goal_int(0),  _goal_int(1));

  std::vector<StatePtr> path = gs.getPath();
  if (path.size() < 1) {
    if(_planner_verbose)
      printf(ANSI_COLOR_RED "Cannot find a path fron [%f, %f, %f] to [%f, %f, %f], Abort!" ANSI_COLOR_RESET "\n",
	  start(0), start(1), start(2), goal(0), goal(1), goal(2));
    _status = -1;
    return false;
  }

  //**** raw path, s --> g
  vec_Vec3f ps;
  for (const auto &it : path)
    ps.push_back(_map_util->intToFloat(Vec3i(it->x, it->y, it->z)));
  _raw_path = ps;
  //std::reverse(std::begin(_raw_path), std::end(_raw_path));

  /*
  if (_goal_outside) // if the second last point is frontier, remove goal
    _raw_path.pop_back();
    */

  _path = removePts(_raw_path);
  _path = crop(_path);
  _path = optimize(_path);
  std::reverse(std::begin(_path), std::end(_path));
  _path = optimize(_path);
  std::reverse(std::begin(_path), std::end(_path));
  return true;
}
