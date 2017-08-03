/**
 * @file nx_jps_3d.h
 * @brief JPS::NXJPS3DUtil Class
 */

#include <jps3d/planner_base.h>
#include <jps3d/nx_planner/jps_3D.h>
#include <iostream>

namespace JPS {
  class NXJPS3DUtil : public PlannerBase {
    public:
      NXJPS3DUtil(bool verbose = false);
      bool plan(const Vec3f &start, const Vec3f &goal, decimal_t eps = 1.0);
    private:
      bool linkToGoal(const Vec3i &pn);

      Vec3i _start_int;
      Vec3i _goal_int;
  };
}
