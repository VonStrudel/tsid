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
#include "tsid/tasks/task-se3-equality.hpp"
#include "tsid/robots/robot-wrapper.hpp"

namespace tsid
{
  namespace tasks
  {
    using namespace math;
    using namespace trajectories;
    using namespace pinocchio;

    TaskTest::TaskTest(const std::string & name,
                                     RobotWrapper & robot,
                                     const std::string & frameName):
      TaskMotion(name, robot),
      m_frame_name(frameName),
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
    }

    void TaskTest::setMask(math::ConstRefVector mask)
    {
      TaskMotion::setMask(mask);
      m_constraint.resize(dim(), (unsigned int)m_J.cols());
    }

    int TaskTest::dim() const
    {
      return (int)m_mask.sum();
    }

    const Vector & TaskTest::Kp() const { return m_Kp; }

    const Vector & TaskTest::Kd() const { return m_Kd; }

    void TaskTest::Kp(ConstRefVector Kp)
    {
      assert(Kp.size()==6);
      m_Kp = Kp;
    }

    void TaskTest::Kd(ConstRefVector Kd)
    {
      assert(Kd.size()==6);
      m_Kd = Kd;
    }

    void TaskTest::setReference(TrajectorySample & ref)
    {
      m_ref = ref;
      vectorToSE3(ref.pos, m_M_ref);
      m_v_ref = Motion(ref.vel);
      m_a_ref = Motion(ref.acc);
    }

    const TrajectorySample & TaskTest::getReference() const
    {
      return m_ref;
    }

    const Vector & TaskTest::position_error() const
    {
      return m_p_error_vec;
    }

    const Vector & TaskTest::velocity_error() const
    {
      return m_v_error_vec;
    }

    const Vector & TaskTest::position() const
    {
      return m_p;
    }

    const Vector & TaskTest::velocity() const
    {
      return m_v;
    }

    const Vector & TaskTest::position_ref() const
    {
      return m_p_ref;
    }

    const Vector & TaskTest::velocity_ref() const
    {
      return m_v_ref_vec;
    }

    const Vector & TaskTest::getDesiredAcceleration() const
    {
      return m_a_des;
    }

    Vector TaskTest::getAcceleration(ConstRefVector dv) const
    {
      return m_constraint.matrix()*dv + m_drift.toVector();
    }

    Index TaskTest::frame_id() const
    {
      return m_frame_id;
    }

    const ConstraintBase & TaskTest::getConstraint() const
    {
      return m_constraint;
    }

    void TaskTest::useLocalFrame(bool local_frame)
    {
      m_local_frame = local_frame;
    }

    const ConstraintBase & TaskTest::compute(const double ,
                                                    ConstRefVector ,
                                                    ConstRefVector ,
                                                    const Data & data)
    {
      SE3 oMi;
      Motion v_frame;
      m_robot.framePosition(data, m_frame_id, oMi);
      m_robot.frameVelocity(data, m_frame_id, v_frame);
      m_robot.frameClassicAcceleration(data, m_frame_id, m_drift);

      // @todo Since Jacobian computation is cheaper in world frame
      // we could do all computations in world frame
      m_robot.frameJacobianLocal(data, m_frame_id, m_J);

      errorInSE3(oMi, m_M_ref, m_p_error);          // pos err in local frame
      m_p_error_vec = m_p_error.toVector();
      SE3ToVector(m_M_ref, m_p_ref);
      SE3ToVector(oMi, m_p);

      // Transformation from local to world
      m_wMl.rotation(oMi.rotation());

      if (m_local_frame) {
        m_v_error = v_frame - m_wMl.actInv(m_v_ref);  // vel err in local frame

        // desired acc in local frame
        m_a_des = - m_Kp.cwiseProduct(m_p_error_vec)
                  - m_Kd.cwiseProduct(m_v_error.toVector())
                  + m_wMl.actInv(m_a_ref).toVector();
      } else {
        m_p_error_vec = m_wMl.toActionMatrix() *   // pos err in local world-oriented frame
            m_p_error.toVector();
        m_v_error = m_wMl.act(v_frame) - m_v_ref;  // vel err in local world-oriented frame

        m_drift = m_wMl.act(m_drift);

        // desired acc in local oriented frame
        m_a_des = - m_Kp.cwiseProduct(m_p_error_vec)
                  - m_Kd.cwiseProduct(m_v_error.toVector())
                  + m_a_ref.toVector();

        // Use an explicit temporary `m_J_rotated` here to avoid allocations.
        m_J_rotated.noalias() = m_wMl.toActionMatrix() * m_J;
        m_J = m_J_rotated;
      }

      m_v_error_vec = m_v_error.toVector();
      m_v_ref_vec = m_v_ref.toVector();
      m_v = v_frame.toVector();


      int idx = 0;
      for (int i = 0; i < 6; i++) {
        if (m_mask(i) != 1.) continue;

        m_constraint.matrix().row(idx) = m_J.row(i);
        m_constraint.vector().row(idx) = (m_a_des - m_drift.toVector()).row(i);

        idx += 1;
      }
      return m_constraint;
    }
  }
}
