//
// Copyright (c) 2017 CNRS, NYU, MPI Tübingen
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

#include "tsid/math/utils.hpp"
#include "tsid/tasks/task-wall.hpp"
#include "tsid/robots/robot-wrapper.hpp"

namespace tsid
{
  namespace tasks
  {
    using namespace math;
    using namespace trajectories;
    using namespace pinocchio;

    TaskWall::TaskWall(const std::string & name,
                                     RobotWrapper & robot,
                                     const std::string & frameName):
      TaskMotion(name, robot),
      m_frame_name(frameName),
      //name rows cols
      m_constraint(name, 6, robot.nv()),
      m_ref(12, 6)
    {
      assert(m_robot.model().existFrame(frameName));
      m_frame_id = m_robot.model().getFrameId(frameName);

      m_v_ref.setZero();
      m_a_ref.setZero();
      m_M_ref.setIdentity();
      m_wMl.setIdentity();
      m_p_error_vec.setZero(6);
      m_v_error_vec.setZero(6);
      m_p.resize(12);
      m_v.resize(6);
      m_p_ref.resize(12);
      m_v_ref_vec.resize(6);
      m_Kp.setZero(6);
      m_Kd.setZero(6);
      m_a_des.setZero(6);
      m_J.setZero(6, robot.nv());
      m_J_rotated.setZero(6, robot.nv());

      m_mask.resize(6);
      m_mask.fill(1.);

      m_local_frame = true;

      m_a = 0.0;
      m_b = 0.0;
      m_c = 0.0;
      m_d = 0.0;
      m_sqrtabc = 0.0;
      m_name = name;
    }

    void TaskWall::setMask(math::ConstRefVector mask)
    {
      TaskMotion::setMask(mask);
      m_constraint.resize(dim(), (unsigned int)m_J.cols());
    }

    int TaskWall::dim() const
    {
      return (int)m_mask.sum();
    }

    const Vector & TaskWall::Kp() const { return m_Kp; }

    const Vector & TaskWall::Kd() const { return m_Kd; }

    void TaskWall::Kp(ConstRefVector Kp)
    {
      assert(Kp.size()==6);
      m_Kp = Kp;
    }

    void TaskWall::Kd(ConstRefVector Kd)
    {
      assert(Kd.size()==6);
      m_Kd = Kd;
    }

    void TaskWall::setReference(TrajectorySample & ref)
    {
      m_ref = ref;
      vectorToSE3(ref.pos, m_M_ref);
      m_v_ref = Motion(ref.vel);
      m_a_ref = Motion(ref.acc);
    }

    const TrajectorySample & TaskWall::getReference() const
    {
      return m_ref;
    }

    const Vector & TaskWall::position_error() const
    {
      return m_p_error_vec;
    }

    const Vector & TaskWall::velocity_error() const
    {
      return m_v_error_vec;
    }

    const Vector & TaskWall::position() const
    {
      return m_p;
    }

    const Vector & TaskWall::velocity() const
    {
      return m_v;
    }

    const Vector & TaskWall::position_ref() const
    {
      return m_p_ref;
    }

    const Vector & TaskWall::velocity_ref() const
    {
      return m_v_ref_vec;
    }

    const Vector & TaskWall::getDesiredAcceleration() const
    {
      return m_a_des;
    }

    Vector TaskWall::getAcceleration(ConstRefVector dv) const
    {
      return m_constraint.matrix()*dv + m_drift.toVector();
    }

    Index TaskWall::frame_id() const
    {
      return m_frame_id;
    }

    const ConstraintBase & TaskWall::getConstraint() const
    {
      return m_constraint;
    }

    void TaskWall::useLocalFrame(bool local_frame)
    {
      m_local_frame = local_frame;
    }

    void TaskWall::set_abcd(double a, double b, double c, double d)
    {
      m_a = a;
      m_b = b;
      m_c = c;
      m_d = d;

      m_sqrtabc = sqrt(a*a+b*b+c*c);
    }

    double TaskWall::l()
    {
      double _l = 0.0;
      double _x = m_p[0];
      double _y = m_p[1];
      double _z = m_p[2];

      _l = m_a*_x + m_b*_y + m_c*_z;
      _l = _l / m_sqrtabc;
      return _l;
    }

    //Jl : 1*nv
    //Changement de distance du mur en fonction des mouvements articulaires
    Vector TaskWall::Jl(ConstRefVector m_J)
    {
      Vector _Jl= Vector(m_robot.nv());
      Vector _n = Vector(m_robot.nv());
      _n(0) = m_a;
      _n(1) = m_b;
      _n(2) = m_c;

      _Jl.setZero(m_robot.nv());
      _Jl =  _n.transpose() * m_J;
      return _Jl;
    }

    const ConstraintBase & TaskWall::compute(const double t,
                                                    ConstRefVector,
                                                    ConstRefVector,
                                                    const Data & data)
    {
      static double old_t = 0.0;
      double dt = t-old_t;
      old_t = t;
      SE3 oMi;

      //Dans robotWrapper.cpp : pinocchio.getJointJacobian
      m_robot.frameJacobianLocal(data, m_frame_id, m_J);

      //Position dans l'espace de la frame qui nous interesse
      m_robot.framePosition(data, m_frame_id, oMi);
      SE3ToVector(oMi, m_p);

      double _l = l();
      static double old_l = 0.0;
      double dl = _l - old_l;
      old_l = _l;

      Vector m_Jl = Jl(m_J);

      Vector lb = Vector(1);
      Vector ub = Vector(1);
                      //rows cols
      Matrix A = Matrix(1, m_robot.nv());//Vecteur ligne

      lb(0) = (2*m_c-_l-dt*dl)/(dt*dt);
      ub(0) = 10000.0;
      //A.row(0) = m_Jl;
      m_constraint = ConstraintInequality(m_name,m_Jl,lb,ub);
      return m_constraint;
    }
  }
}
