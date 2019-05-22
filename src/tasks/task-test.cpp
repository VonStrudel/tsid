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
      //                rows     cols
      m_constraint(name, 1, robot.nv()),
      //             robot.nv : dimension of the velocity vector space
      m_ref(12, 6)
    {
      assert(m_robot.model().existFrame(frameName));

      //numéro de la frame sur laquelle on applique la tache
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


    //notes :
    //Dans fwd.hpp (inclu dans constraint-x.hpp inclu dans taskTest.hpp)
    //typedef Eigen::Matrix<Scalar,Eigen::Dynamic,1> Vector;
    //typedef Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> Matrix;


    void TaskTest::setABCD(int a, int b, int c, int d)
    {
      m_a = a;
      m_b = b;
      m_c = c;
      m_d = d;
      m_sqrtabc = sqrt(a*a+b*b+c*c);
      m_n = Vector3d(a,b,c);
    }

    /*
    static pinocchio::SE3 position(const Robot & self, const pinocchio::Data & data, const pinocchio::Model::JointIndex & index){
      return self.position(data, index);
    }

    //in utils.hpp
    //Convert the input SE3 object to a 7D vector of floats [X,Y,Z,Q1,Q2,Q3,Q4].
    void SE3ToXYZQUAT(const pinocchio::SE3 & M, RefVector xyzQuat);
    //RefVector is defined in fwd.hpp as
    typedef Eigen::Ref<Vector>              RefVector;
    */
    inline double TaskTest::getL(const Data & data)
    //L : wall distance relative to hand
    {
      SE3 frame_position;
      Vector xyzquat_frame_position(7);

      frame_position = m_robot.position(data, m_frame_id); //coordonnées x,y,z de la frame
      SE3ToXYZQUAT(frame_position, xyzquat_frame_position);

      return (mxyzquat_frame_position_x[0]*a
        +xyzquat_frame_position[1]*b
        +xyzquat_frame_position[2]*c+d)
        /m_sqrtabc;
    }



    const ConstraintBase & TaskTest::compute(const double t,
                                                    ConstRefVector q ,
                                                    ConstRefVector v,
                                                    const Data & data)
    {

      ERREUR FLAGRANTE
      static double prev_t=0;
      static double prev_l=0;

      //l : distance between hand and wall
      //\dot l : speed of hand relative to wall in normal vector direction
      double l, dl;
      Vector lb(1);//define size 1,1
      Vector ub(1);
      //Jl : jacobian of wall distance relative to q
      Matrix Jl;//define size 1 row, nv cols
      Matrix A;//define size 1 row, nv cols

      m_robot.frameJacobianLocal(data, m_frame_id, m_J);

      //Constraint structure
      //lb <= A <= ub
      Jl = m_n.transpose()*m_J;// getJl(q);
      l = getL(q);
      dl = l-prev_l;
      dt = t-prev_t;
      A = Jl;
      lb(0,0) = (2*m_c-l-dt*dl)/(dt*dt);//-dJl*v;
      ub(0,0) =
      m_constraint.setMatrix(Jl);
      m_constraint.setUpperBound(ub);
      m_constraint.setLowerBound(lb);
      prev_t = t;
      prev_l = l;

      return m_constraint;

    }
  }
}
