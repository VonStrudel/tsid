//
// Copyright (c) 2018 CNRS
//
// This file is part of tsid
// tsid is free software: you can redistribute it
// and/or modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation, either version
// 3 of the License, or (at your option) any later version.
// tsid is distributed in the hope that it will be
// useful, but WITHOUT ANY WARRANTY; without even the implied warranty
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// General Lesser Public License for more details. You should have
// received a copy of the GNU Lesser General Public License along with
// tsid If not, see
// <http://www.gnu.org/licenses/>.
//

#ifndef __tsid_python_task_wall_hpp__
#define __tsid_python_task_wall_hpp__

#include <boost/python.hpp>
#include <string>
#include <eigenpy/eigenpy.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "tsid/tasks/task-wall.hpp"
#include "tsid/robots/robot-wrapper.hpp"
#include "tsid/trajectories/trajectory-base.hpp"
#include "tsid/math/constraint-inequality.hpp"
#include "tsid/math/constraint-base.hpp"
namespace tsid
{
  namespace python
  {
    namespace bp = boost::python;

    template<typename TaskWall>
    struct TaskWallPythonVisitor
    : public boost::python::def_visitor< TaskWallPythonVisitor<TaskWall> >
    {

      template<class PyClass>


      void visit(PyClass& cl) const
      {
        cl
        .def(bp::init<std::string, robots::RobotWrapper &, std::string> ((bp::arg("name"), bp::arg("robot"), bp::arg("framename")), "Default Constructor"))
        .add_property("dim", &TaskWall::dim, "return dimension size")
        .def("setReference", &TaskWallPythonVisitor::setReference, bp::arg("ref"))
        .add_property("getDesiredAcceleration", bp::make_function(&TaskWallPythonVisitor::getDesiredAcceleration, bp::return_value_policy<bp::copy_const_reference>()), "Return Acc_desired")
        .def("getAcceleration", &TaskWallPythonVisitor::getAcceleration, bp::arg("dv"))
        .add_property("position_error", bp::make_function(&TaskWallPythonVisitor::position_error, bp::return_value_policy<bp::copy_const_reference>()))
        .add_property("velocity_error", bp::make_function(&TaskWallPythonVisitor::velocity_error, bp::return_value_policy<bp::copy_const_reference>()))
        .add_property("position", bp::make_function(&TaskWallPythonVisitor::position, bp::return_value_policy<bp::copy_const_reference>()))
        .add_property("velocity", bp::make_function(&TaskWallPythonVisitor::velocity, bp::return_value_policy<bp::copy_const_reference>()))
        .add_property("position_ref", bp::make_function(&TaskWallPythonVisitor::position_ref, bp::return_value_policy<bp::copy_const_reference>()))
        .add_property("velocity_ref", bp::make_function(&TaskWallPythonVisitor::velocity_ref, bp::return_value_policy<bp::copy_const_reference>()))
        .add_property("Kp", bp::make_function(&TaskWallPythonVisitor::Kp, bp::return_value_policy<bp::copy_const_reference>()))
        .add_property("Kd", bp::make_function(&TaskWallPythonVisitor::Kd, bp::return_value_policy<bp::copy_const_reference>()))
        .def("setKp", &TaskWallPythonVisitor::setKp, bp::arg("Kp"))
        .def("setKd", &TaskWallPythonVisitor::setKd, bp::arg("Kd"))
        .def("useLocalFrame", &TaskWallPythonVisitor::useLocalFrame, bp::arg("local_frame"))
        .def("setMask", &TaskWallPythonVisitor::setMask, bp::arg("mask"))
        .def("compute", &TaskWallPythonVisitor::compute, bp::args("t", "q", "v", "data"))
        .def("getConstraint",  &TaskWallPythonVisitor::getConstraint)
        .def("set_abcd", &TaskWallPythonVisitor::set_abcd, bp::args("a","b","c","d"))
        .add_property("frame_id", &TaskWall::frame_id, "frame id return")
        .add_property("name", &TaskWallPythonVisitor::name)
        ;
      }
       static std::string name(TaskWall & self){
        std::string name = self.name();
        return name;
      }
      static math::ConstraintInequality compute(TaskWall & self, const double t, const Eigen::VectorXd & q, const Eigen::VectorXd & v, const pinocchio::Data & data){
        self.compute(t, q, v, data);
        math::ConstraintInequality cons(self.getConstraint().name(), self.getConstraint().matrix(), self.getConstraint().lowerBound(),self.getConstraint().upperBound());
        return cons;
      }
      static math::ConstraintInequality getConstraint(const TaskWall & self){
        math::ConstraintInequality cons(self.getConstraint().name(), self.getConstraint().matrix(), self.getConstraint().lowerBound(),self.getConstraint().upperBound());
        return cons;
      }
      static void setReference(TaskWall & self, trajectories::TrajectorySample & ref){
        self.setReference(ref);
      }
      static const Eigen::VectorXd & getDesiredAcceleration(const TaskWall & self){
        return self.getDesiredAcceleration();
      }
      static Eigen::VectorXd getAcceleration (TaskWall & self, const Eigen::VectorXd dv){
        return self.getAcceleration(dv);
      }
      static const Eigen::VectorXd & position_error(const TaskWall & self){
        return self.position_error();
      }
      static const Eigen::VectorXd & velocity_error(const TaskWall & self){
        return self.velocity_error();
      }
      static const Eigen::VectorXd & position (const TaskWall & self){
        return self.position();
      }
      static const Eigen::VectorXd & velocity (const TaskWall & self){
        return self.velocity();
      }
      static const Eigen::VectorXd & position_ref (const TaskWall & self){
        return self.position_ref();
      }
      static const Eigen::VectorXd & velocity_ref (const TaskWall & self){
        return self.velocity_ref();
      }
      static const Eigen::VectorXd & Kp (TaskWall & self){
        return self.Kp();
      }
      static const Eigen::VectorXd & Kd (TaskWall & self){
        return self.Kd();
      }
      static void setKp (TaskWall & self, const::Eigen::VectorXd Kp){
        return self.Kp(Kp);
      }
      static void setKd (TaskWall & self, const::Eigen::VectorXd Kv){
        return self.Kd(Kv);
      }
      static void useLocalFrame (TaskWall & self, const bool local_frame) {
        self.useLocalFrame(local_frame);
      }
      static void setMask (TaskWall & self, const::Eigen::VectorXd mask) {
        self.setMask(mask);
      }
      static Eigen::VectorXd frame_id (TaskWall & self){
        return self.frame_id();
      }
      static void set_abcd(TaskWall & self, double a, double b, double c, double d)
      {
        self.set_abcd(a,b,c,d);
      }
      static void expose(const std::string & class_name)
      {
        std::string doc = "TaskWall info.";
        bp::class_<TaskWall>(class_name.c_str(),
                          doc.c_str(),
                          bp::no_init)
        .def(TaskWallPythonVisitor<TaskWall>());

      //  bp::register_ptr_to_python< boost::shared_ptr<math::ConstraintBase> >();
      }
    };
  }
}


#endif // ifndef __tsid_python_task_wall_hpp__
