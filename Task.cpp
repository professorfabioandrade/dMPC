//***************************************************************************
// Copyright 2007-2017 Universidade do Porto - Faculdade de Engenharia      *
// Laboratório de Sistemas e Tecnologia Subaquática (LSTS)                  *
//***************************************************************************
// This file is part of DUNE: Unified Navigation Environment.               *
//                                                                          *
// Commercial Licence Usage                                                 *
// Licencees holding valid commercial DUNE licences may use this file in    *
// accordance with the commercial licence agreement provided with the       *
// Software or, alternatively, in accordance with the terms contained in a  *
// written agreement between you and Faculdade de Engenharia da             *
// Universidade do Porto. For licensing terms, conditions, and further      *
// information contact lsts@fe.up.pt.                                       *
//                                                                          *
// Modified European Union Public Licence - EUPL v.1.1 Usage                *
// Alternatively, this file may be used under the terms of the Modified     *
// EUPL, Version 1.1 only (the "Licence"), appearing in the file LICENCE.md *
// included in the packaging of this file. You may not use this work        *
// except in compliance with the Licence. Unless required by applicable     *
// law or agreed to in writing, software distributed under the Licence is   *
// distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF     *
// ANY KIND, either express or implied. See the Licence for the specific    *
// language governing permissions and limitations at                        *
// https://github.com/LSTS/dune/blob/master/LICENCE.md and                  *
// http://ec.europa.eu/idabc/eupl.html.                                     *
//***************************************************************************
// Author: FabioAndrade                                                     *
//***************************************************************************

// DUNE headers.
#include <DUNE/DUNE.hpp>
#include "kernel.h"

namespace MultiAgent
{
  namespace UAVagent
  {
    using DUNE_NAMESPACES;

    struct Task: public DUNE::Tasks::Task
    {
      uint8_t index;

      bool controls_received[I] = { 0 };        //boolean array to monitor if the agent got the control inputs of other agents
      bool all_controls_received;               //a boolean variable to start the MPC-PSO when all control inputs of other agents are received

      struct agent agents[I];                   //struct with state and control inputs of agents
      struct point2d r[Nb_x][Nb_y];             //x,y position of grid 
      
      //Create IMC messages
      IMC::DesiredSpeed dspeed;                 //Create IMC message that sends the desired speed
      IMC::DesiredRoll droll;                   //Create IMC message that sends the desired roll
      IMC::ControlLoops cloops;                 //Create IMC message that sends the desired control loops

      uint64_t current_time;                    //Time right after sending the control inputs to the autopilot
      uint64_t last_time;                       //Last time that the control inputs were sent to the autopilot
      uint64_t before_PSO;                      //Time before the PSO calculation
      uint64_t after_PSO;                       //Time after the PSO calculation
      uint64_t state_time[I];                   //Time where the current states were gotten (for forward euler method)
      double Dt_tot = Dt_PSO;                   //Delta time between the last time that the control inputs were sent and the current
      double Dt_PSO = 0.25;                     //Delta time between the before and after PSO
      double Dt_state[I] = { 0 };              //Delta time between the state time  of the agents and the time before PSO

      bool first = 1;                           //control boolean to start and stop the MPC-PSO

      //! Constructor.
      //! @param[in] name task name.
      //! @param[in] ctx context.
      Task(const std::string& name, Tasks::Context& ctx):
        DUNE::Tasks::Task(name, ctx)
      {
        param("Index", index);                    //Agent param (changable in the .ini file)

        bind<IMC::Target>(this);                  //To bind to the Target IMC message (received from the AUV)
        bind<IMC::EstimatedState>(this);          //To bind the EstimatedState IMC message (that has the current state)
        bind<IMC::IndicatedSpeed>(this);          //To bind the IndicatedSpeed IMC message (that has the airspeed)
        bind<IMC::multiagent>(this);
      }

      void
      consume(const IMC::IndicatedSpeed* airspeed)
      {
        agents[index].v_a0 = airspeed->value;   //records the airspeed value in the agent current state vector
      }

      void
      consume(const IMC::EstimatedState* m_estimatedstate) //gets indicated airspeed to use as current airspeed
      {

        //Reference-LLH of agent (it is not its llh, but the reference llh used by the agent)
        agents[index].lat = m_estimatedstate->lat;        
        agents[index].lon = m_estimatedstate->lon;
        agents[index].height = m_estimatedstate->height;
        //Convert local NED (that uses the Reference-LLH) to Agent-LLH (agent's llh)
        WGS84::displace(m_estimatedstate->x,m_estimatedstate->y,m_estimatedstate->z,&agents[index].lat,&agents[index].lon,&agents[index].height);
        //Convert from Agent-LLH to global NED (which uses the MPC-PSO LLH reference)
        WGS84::displacement(rlat,rlon,0,agents[index].lat,agents[index].lon,agents[index].height,&agents[index].x[0],&agents[index].y[0],&agents[index].z);
        agents[index].psi[0] = m_estimatedstate->psi; //current yaw angle
        agents[index].phi0 = m_estimatedstate->phi; //current roll angle
        state_time[index] = Clock::getMsec();
        //float speed_g = Math::norm(Math::norm(m_estimatedstate->vx,m_estimatedstate->vy),m_estimatedstate->vz);//ts.speed;
        //float speed_g2 = sqrt( pow(agents[index].v_a0*cos(agents[index].psi[0]) + agents[index].v_w*cos(agents[index].psi_w), 2) + pow(agents[index].v_a0*sin(agents[index].psi[0]) + agents[index].v_w*sin(agents[index].psi_w),2) );

        inf("\t %5.2f\t %5.2f\t %5.2f\t %5.2f\t %5.2f\t %5.2f\t %5.2f\t %5.2f\t",
          agents[index].x[0],agents[index].y[0],agents[index].psi[0]*180/pi,agents[index].v_a0,agents[index].phi0*180/pi,sqrt(pow(agents[index].x[0] - agents[0].x[0], 2) + pow(agents[index].y[0] - agents[0].y[0], 2)),sqrt(pow(agents[index].x[0] - agents[1].x[0], 2) + pow(agents[index].y[0] - agents[1].y[0], 2)),sqrt(pow(agents[index].x[0] - agents[2].x[0], 2) + pow(agents[index].y[0] - agents[2].y[0], 2)));
      
        //check if Nb_x x Nb_y cells are visited
        int x_idx = (int)agents[index].x[0]/res;
        int y_idx = (int)agents[index].y[0]/res;
        //debug("F = %f",calculate_F(agents[index].x[0],agents[index].y[0],agents[index].B,r,phi_out));
        //debug("Cell ----- [%d,%d]",x_idx,y_idx);
        for (int i = (x_idx - R/res); i <= (x_idx + R/res); i++)
        {
          for (int j = (y_idx - R/res); j <= (y_idx + R/res); j++)
          {
            //debug("i=%d,j%d \n",i,j);
            if ( (i>=0) && (i<Nb_x) && (j>=0) && (j<Nb_y) && check_if_grid_is_inside_R(agents[index].x[0],agents[index].y[0],r_out[i][j]))
            {
              //debug("i,j[%d,%d] idx [%d,%d] B %d",i,j,x_idx,y_idx,agents[index].B[i][j]);
              if (agents[index].B[i][j] == 0) debug("----------- VISITED [%d,%d] phi=%f [%d,%d] ---------------",i,j,phi_out[i][j],x_idx,y_idx);
              //if (agents[index].B[i][j] == 1) debug(" ALREADY ---------- [%d,%d] phi=%f [%d,%d] ---------------",i,j,phi_out[i][j],x_idx,y_idx);
              
              agents[index].B[i][j] = 1;
            }
          }
        }

      }

      void
      consume(const IMC::Target* trg)
      {
        if ((trg->z == 1) && (first == 1)) //start the MPC-PSO and enable speed and roll control and send loiter state
        {
          cloops.enable = IMC::ControlLoops::CL_ENABLE;
          cloops.mask = IMC::CL_SPEED;
          dispatch(cloops);
          cloops.mask = IMC::CL_ROLL;
          dispatch(cloops);
          first = 0;
          inf("Start MPC-PSO");
          //To tell the Ardupilot to keep the current loiter speed and roll angle
          dspeed.value = agents[index].v_a0;
          droll.value = agents[index].phi0;
          dispatch(dspeed);
          dispatch(droll);
          //Send to other agents
          for (int k = 0; k < N; k++) agents[index].u_v[k] = agents[index].v_a0;
          for (int k = 0; k < N; k++) agents[index].u_phi[k] = agents[index].phi0;
          send_self(index, agents);

          debug("First controls dispatched (Loiter) at %ld",Clock::getMsec());

          //Load the probabilistic map of the grids
        }

        if ((trg->z == 0) && (first == 0)) //disable MPC-PSO
        {
          cloops.enable = IMC::ControlLoops::CL_DISABLE;
          cloops.mask = IMC::CL_SPEED;
          dispatch(cloops);
          cloops.mask = IMC::CL_ROLL;
          dispatch(cloops);
          first = 1;
          inf("Stop MPC-PSO");
          IMC::IdleManeuver Idle;   //creates idle maneuver IMC message
          dispatch(Idle);           //dispatches idle maneuver
        }

        if ((trg->z == 0) && (first == 1)) {} //do nothing
        if ((trg->z == 1) && (first == 0)) {} //do nothing
      }

      //! Update internal state with new parameter values.
      void
      onUpdateParameters(void)
      {
      }

      //! Reserve entity identifiers.
      void
      onEntityReservation(void)
      {
      }

      //! Resolve entity names.
      void
      onEntityResolution(void)
      {
      }

      //! Acquire resources.
      void
      onResourceAcquisition(void)
      {
      }

      //! Initialize resources.
      void
      onResourceInitialization(void)
      {
        loadphi();

        for (int i = 0; i < Nb_x; i++)
        {
          for (int j = 0; j < Nb_y; j++)
          {
            r_out[i][j].x = (x_min + res/2) + res*i;
            r_out[i][j].y = (y_min + res/2) + res*j;
          }
        } 
      }

      //! Release resources.
      void
      onResourceRelease(void)
      {
      }

      //! Main loop.
      void
      onMain(void)
      {
        while (!stopping())
        {
          if (first == 0)
          {
            controls_received[index] = 1; //To say that the current agent is ready
            all_controls_received = 1;    //To restart the boolean
            for (int i = 0; i < I; i++) all_controls_received = all_controls_received && controls_received[i]; //Check if the control inputs of all agents were received
            if (all_controls_received) 
            {
              before_PSO = Clock::getMsec();
              //debug("Before up \t %5.2f\t %5.2f\t %5.2f\t %5.2f\t %5.2f\t",
              //agents[index].x[0],agents[index].y[0],agents[index].psi[0],agents[index].v_a0,agents[index].phi0*180/pi);
              for (int i = 0; i < I; i++) Dt_state[i] = ((double)before_PSO - (double)state_time[i])/1000.0;
              //debug("dt0=%lf,dt1=%lf,dt2=%lf",Dt_state[0],Dt_state[1],Dt_state[2]);
              for (int i = 0; i < I; i++) update_initialstate_with_delay(i,agents,(Dt_PSO + Dt_state[i]));
              for (int i = 0; i < I; i++) forward_euler(agents[i].x, agents[i].y, agents[i].psi, agents[i].v_w, agents[i].psi_w, agents[i].u_phi, agents[i].u_v, Dt_MPC);
              //debug("After up  \t %5.2f\t %5.2f\t %5.2f\t %5.2f\t %5.2f\t",
              //agents[index].x[0],agents[index].y[0],agents[index].psi[0],agents[index].v_a0,agents[index].phi0*180/pi);

              PSO(index, agents); //Calls the MPC-PSO algorithm, which also updates the states of the agents
              after_PSO = Clock::getMsec();
              Dt_PSO = ((double)after_PSO - (double)before_PSO)/1000.0;
              //v_a
              dspeed.value = agents[index].u_v[0];
              //phi
              droll.value = agents[index].u_phi[0];
              //constrains
              if (dspeed.value < v_min) dspeed.value = v_min; //(Eq. 05)
              if (dspeed.value > v_max) dspeed.value = v_max;
              if (droll.value < phi_min) droll.value = phi_min; //(Eq. 06)
              if (droll.value > phi_max) droll.value = phi_max;
              //Dispatch control commands
              dispatch(dspeed);
              dispatch(droll);
              current_time = Clock::getMsec();
              Dt_tot = ((double)current_time - (double)last_time)/1000.0;
              last_time = current_time;
              //debug("Dt_tot = %lf , Dt_PSO = %lf",Dt_tot,Dt_PSO);
              inf("Dt_PSO=%f,phi=%f->%f, v_a=%f->%f",Dt_PSO,agents[index].phi0*180/pi,droll.value*180/pi,agents[index].v_a0,dspeed.value);

              


                ////////////////////////////////////////////////////////////
                  //update the cells that the agents are planning to visit

              float controls[N*2];
              for (int k=0;k<N;k++) controls[k] = agents[index].u_phi[k];
              for (int k=0;k<N;k++) controls[k+N] = agents[index].u_v[k];


                  for (int i = 0; i < I; i++)
                  {
                    if (1)
                    {
                      for (int k = 0; k < (N+1); k++)
                      { 
                        //check if cells are going to be visited
                        int x_idx = (int)agents[i].x[k]/res;
                        int y_idx = (int)agents[i].y[k]/res;
                        for (int l = (x_idx - R/res); l <= (x_idx + R/res); l++)
                        {
                          for (int j = (y_idx - R/res); j <= (y_idx + R/res); j++)
                          {
                            //printf("x=%f,y=%f \n",x1[k],y1[k]);
                            if ( (l>=0) && (l<Nb_x) && (j>=0) && (j<Nb_y) && check_if_grid_is_inside_R(agents[i].x[k],agents[i].y[k],r[l][j]) )
                            {
                              agents[i].b[l][j] = 1;
                            }
                          }
                        }
                      }
                    }
                  }

              debug("cost=%f",cost_function(index,controls,agents));


              //zero all the cells that the agents are planning to visit
              for (int i = 0; i < I; i++)
              {
                for (int l = 0; l < Nb_x; l++)
                {
                  for (int j = 0; j < Nb_y; j++)
                  {
                    agents[i].b[l][j] = 0;
                  }
                }
              }
              //////////////////////////////////////////////////////

              send_self(index, agents);
              trace("State and Control Inputs dispatched at %ld", Clock::getMsec());
              //To zero the boolean
              for (int i = 0; i < I; i++) controls_received[i] = 0;
            }
            else 
            {
              trace("All controls not received yet at %ld",Clock::getMsec());
            }
          }
        waitForMessages(1.0);
        }
      }

      void update_initialstate_with_delay(int index_, struct agent *agents_, float Dt)
      {
        float v_g = sqrt( pow((agents_[index_].v_a0 + agents_[index_].u_v[0])/2*cos(agents_[index_].psi[0]) + agents_[index_].v_w*cos(agents_[index_].psi_w), 2) + pow((agents_[index_].v_a0 + agents_[index_].u_v[0])/2*sin(agents_[index_].psi[0]) + agents_[index_].v_w*sin(agents_[index_].psi_w),2) );
        float chi = atan2( ( (agents_[index_].v_a0 + agents_[index_].u_v[0])/2*sin(agents_[index_].psi[0]) + agents_[index_].v_w*sin(agents_[index_].psi_w) ) , ( (agents_[index_].v_a0 + agents_[index_].u_v[0])/2*cos(agents_[index_].psi[0]) + agents_[index_].v_w*cos(agents_[index_].psi_w) ));
        if (chi < 0) chi = chi + 2*pi;
        //x
        float dx = v_g*cos(chi);
        agents_[index_].x[0] = agents_[index_].x[0] + dx*Dt;
        //y
        float dy = v_g*sin(chi);
        agents_[index_].y[0] = agents_[index_].y[0] + dy*Dt;
        //chi
        float dchi = g*tan((agents_[index_].phi0 + agents_[index_].u_phi[0])/2)*cos(chi - agents_[index_].psi[0])/v_g;
        chi = chi + dchi*Dt;
        if (chi > 2*pi) chi = chi - 2*pi;
        //psi
        agents_[index_].psi[0] = chi - asin(agents_[index_].v_w/agents_[index_].v_a0*sin(agents_[index_].psi_w - chi));
        if (agents_[index_].psi[0] > pi) agents_[index_].psi[0] = agents_[index_].psi[0] -2*pi;

        //For broadcast information, it is needed to reduce the first planned controls, discounting the Dt from the Dt_MPC
        agents_[index_].u_v[0] = agents_[index_].u_v[0]*(Dt_MPC - Dt);
      }

      void send_self(int index_, struct agent agents_[])
      {
        IMC::multiagent magt_self;

        magt_self.lat_0 = agents_[index_].lat;
        magt_self.lon_0 = agents_[index_].lon;
        magt_self.height_0 = agents_[index_].height;

        magt_self.psi_0 = agents_[index_].psi[0];
        magt_self.v_a_0 = agents_[index_].v_a0;
        magt_self.phi_0 = agents_[index_].phi0;

        magt_self.u_v_01 = agents_[index_].u_v[0];
        magt_self.u_v_02 = agents_[index_].u_v[1];
        magt_self.u_v_03 = agents_[index_].u_v[2];
        magt_self.u_v_04 = agents_[index_].u_v[3];
        magt_self.u_v_05 = agents_[index_].u_v[4];
        magt_self.u_v_06 = agents_[index_].u_v[5];
        magt_self.u_v_07 = agents_[index_].u_v[6];
        magt_self.u_v_08 = agents_[index_].u_v[7];
        magt_self.u_v_09 = agents_[index_].u_v[8];
        magt_self.u_v_10 = agents_[index_].u_v[9];
        magt_self.u_v_11 = agents_[index_].u_v[10];
        magt_self.u_v_12 = agents_[index_].u_v[11];
        magt_self.u_v_13 = agents_[index_].u_v[12];
        magt_self.u_v_14 = agents_[index_].u_v[13];
        magt_self.u_v_15 = agents_[index_].u_v[14];
        magt_self.u_v_16 = agents_[index_].u_v[15];
        magt_self.u_v_17 = agents_[index_].u_v[16];
        magt_self.u_v_18 = agents_[index_].u_v[17];
        magt_self.u_v_19 = agents_[index_].u_v[18];
        magt_self.u_v_20 = agents_[index_].u_v[19];
        magt_self.u_v_21 = agents_[index_].u_v[20];
        magt_self.u_v_22 = agents_[index_].u_v[21];
        magt_self.u_v_23 = agents_[index_].u_v[22];
        magt_self.u_v_24 = agents_[index_].u_v[23];
        magt_self.u_v_25 = agents_[index_].u_v[24];
        magt_self.u_v_26 = agents_[index_].u_v[25];
        magt_self.u_v_27 = agents_[index_].u_v[26];
        magt_self.u_v_28 = agents_[index_].u_v[27];
        magt_self.u_v_29 = agents_[index_].u_v[28];
        magt_self.u_v_30 = agents_[index_].u_v[29];

        magt_self.u_phi_01 = agents_[index_].u_phi[0];
        magt_self.u_phi_02 = agents_[index_].u_phi[1];
        magt_self.u_phi_03 = agents_[index_].u_phi[2];
        magt_self.u_phi_04 = agents_[index_].u_phi[3];
        magt_self.u_phi_05 = agents_[index_].u_phi[4];
        magt_self.u_phi_06 = agents_[index_].u_phi[5];
        magt_self.u_phi_07 = agents_[index_].u_phi[6];
        magt_self.u_phi_08 = agents_[index_].u_phi[7];
        magt_self.u_phi_09 = agents_[index_].u_phi[8];
        magt_self.u_phi_10 = agents_[index_].u_phi[9];
        magt_self.u_phi_11 = agents_[index_].u_phi[10];
        magt_self.u_phi_12 = agents_[index_].u_phi[11];
        magt_self.u_phi_13 = agents_[index_].u_phi[12];
        magt_self.u_phi_14 = agents_[index_].u_phi[13];
        magt_self.u_phi_15 = agents_[index_].u_phi[14];
        magt_self.u_phi_16 = agents_[index_].u_phi[15];
        magt_self.u_phi_17 = agents_[index_].u_phi[16];
        magt_self.u_phi_18 = agents_[index_].u_phi[17];
        magt_self.u_phi_19 = agents_[index_].u_phi[18];
        magt_self.u_phi_20 = agents_[index_].u_phi[19];
        magt_self.u_phi_21 = agents_[index_].u_phi[20];
        magt_self.u_phi_22 = agents_[index_].u_phi[21];
        magt_self.u_phi_23 = agents_[index_].u_phi[22];
        magt_self.u_phi_24 = agents_[index_].u_phi[23];
        magt_self.u_phi_25 = agents_[index_].u_phi[24];
        magt_self.u_phi_26 = agents_[index_].u_phi[25];
        magt_self.u_phi_27 = agents_[index_].u_phi[26];
        magt_self.u_phi_28 = agents_[index_].u_phi[27];
        magt_self.u_phi_29 = agents_[index_].u_phi[28];
        magt_self.u_phi_30 = agents_[index_].u_phi[29];

        magt_self.i = index_;  //write agent index in the message

        magt_self.data.assign((char*)agents_[index_].B,(char*)agents_[index_].B + sizeof(agents_[index_].B));

        dispatch(magt_self);  //dispatch message
      }

      void
      consume(const IMC::multiagent* magt)
      {
        //debug("Message received from agent %i",magt->i);
        //Agent-LLH of agent i
        agents[magt->i].lat = magt->lat_0;
        agents[magt->i].lon = magt->lon_0;
        agents[magt->i].height = magt->height_0;
        //trace("AGENT %d REF LLH: lat=%lf,lon=%lf,height=%f",magt->i,agents[magt->i].lat,agents[magt->i].lon,agents[magt->i].height);
        //Convert from Agent-LLH to MPC-PSO GLOBAL NED Frame
        WGS84::displacement(rlat,rlon,0.0,agents[magt->i].lat,agents[magt->i].lon,agents[magt->i].height,&agents[magt->i].x[0],&agents[magt->i].y[0],&agents[magt->i].z);
        //trace("AGENT %d GLOBAL NED: x=%f,y=%f,h=%f",magt->i,agents[magt->i].x[0],agents[magt->i].y[0],agents[magt->i].z);
        agents[magt->i].psi[0] = magt->psi_0;
        agents[magt->i].v_a0 = magt->v_a_0;
        agents[magt->i].phi0 = magt->phi_0;
        //trace("AGENT %d: psi=%lf,v_a=%f,phi=%f",magt->i,agents[magt->i].psi[0],agents[magt->i].v_a[0],agents[magt->i].phi[0]);
        agents[magt->i].u_v[0] = magt->u_v_01;
        agents[magt->i].u_v[1] = magt->u_v_02;
        agents[magt->i].u_v[2] = magt->u_v_03;
        agents[magt->i].u_v[3] = magt->u_v_04;
        agents[magt->i].u_v[4] = magt->u_v_05;
        agents[magt->i].u_v[5] = magt->u_v_06;
        agents[magt->i].u_v[6] = magt->u_v_07;
        agents[magt->i].u_v[7] = magt->u_v_08;
        agents[magt->i].u_v[8] = magt->u_v_09;
        agents[magt->i].u_v[9] = magt->u_v_10;
        agents[magt->i].u_v[10] = magt->u_v_11;
        agents[magt->i].u_v[11] = magt->u_v_12;
        agents[magt->i].u_v[12] = magt->u_v_13;
        agents[magt->i].u_v[13] = magt->u_v_14;
        agents[magt->i].u_v[14] = magt->u_v_15;
        agents[magt->i].u_v[15] = magt->u_v_16;
        agents[magt->i].u_v[16] = magt->u_v_17;
        agents[magt->i].u_v[17] = magt->u_v_18;
        agents[magt->i].u_v[18] = magt->u_v_19;
        agents[magt->i].u_v[19] = magt->u_v_20;
        agents[magt->i].u_v[20] = magt->u_v_21;
        agents[magt->i].u_v[21] = magt->u_v_22;
        agents[magt->i].u_v[22] = magt->u_v_23;
        agents[magt->i].u_v[23] = magt->u_v_24;
        agents[magt->i].u_v[24] = magt->u_v_25;
        agents[magt->i].u_v[25] = magt->u_v_26;
        agents[magt->i].u_v[26] = magt->u_v_27;
        agents[magt->i].u_v[27] = magt->u_v_28;
        agents[magt->i].u_v[28] = magt->u_v_29;
        agents[magt->i].u_v[29] = magt->u_v_30;

        agents[magt->i].u_phi[0] = magt->u_phi_01;
        agents[magt->i].u_phi[1] = magt->u_phi_02;
        agents[magt->i].u_phi[2] = magt->u_phi_03;
        agents[magt->i].u_phi[3] = magt->u_phi_04;
        agents[magt->i].u_phi[4] = magt->u_phi_05;
        agents[magt->i].u_phi[5] = magt->u_phi_06;
        agents[magt->i].u_phi[6] = magt->u_phi_07;
        agents[magt->i].u_phi[7] = magt->u_phi_08;
        agents[magt->i].u_phi[8] = magt->u_phi_09;
        agents[magt->i].u_phi[9] = magt->u_phi_10;
        agents[magt->i].u_phi[10] = magt->u_phi_11;
        agents[magt->i].u_phi[11] = magt->u_phi_12;
        agents[magt->i].u_phi[12] = magt->u_phi_13;
        agents[magt->i].u_phi[13] = magt->u_phi_14;
        agents[magt->i].u_phi[14] = magt->u_phi_15;
        agents[magt->i].u_phi[15] = magt->u_phi_16;
        agents[magt->i].u_phi[16] = magt->u_phi_17;
        agents[magt->i].u_phi[17] = magt->u_phi_18;
        agents[magt->i].u_phi[18] = magt->u_phi_19;
        agents[magt->i].u_phi[19] = magt->u_phi_20;
        agents[magt->i].u_phi[20] = magt->u_phi_21;
        agents[magt->i].u_phi[21] = magt->u_phi_22;
        agents[magt->i].u_phi[22] = magt->u_phi_23;
        agents[magt->i].u_phi[23] = magt->u_phi_24;
        agents[magt->i].u_phi[24] = magt->u_phi_25;
        agents[magt->i].u_phi[25] = magt->u_phi_26;
        agents[magt->i].u_phi[26] = magt->u_phi_27;
        agents[magt->i].u_phi[27] = magt->u_phi_28;
        agents[magt->i].u_phi[28] = magt->u_phi_29;
        agents[magt->i].u_phi[29] = magt->u_phi_30;

        int k = 0;
        for (int i = 0; i < Nb_x; i++)
        {
          for (int j = 0; j < Nb_x; j++)
          {
            agents[magt->i].B[i][j] = magt->data[k];
            k = k + 1;
          }
        }

        controls_received[magt->i] = 1;
        state_time[magt->i] = Clock::getMsec();
        trace("Controls received from agent %d at %ld",magt->i,state_time[magt->i]);
      }

      void loadphi() {
        int row, col;
        const int rows = Nb_x;
        const int cols = Nb_y;
        std::ifstream in("/home/fabio/uavlab/dune/src/MultiAgent/UAVagent/phi.txt");
        
        if (in.is_open()) {
          std::cout << "File phi opened";
        }
        else {
          std::cout << "Error opening file phi";
        }

        for (row = 0; row < rows; row++) {
          for (col = 0; col < cols; col++) {
            in >> phi_out[row][col];
            //printf("%f\n",y[row][col]);
          }
        }
        in.close();
      }

    };
  }
}

DUNE_TASK

float cost_function(int index, float controls[], struct agent agents[])
{
    for (int k = 0; k < N; k++) agents[index].u_phi[k] = controls[k];
    for (int k = 0; k < N; k++) agents[index].u_v[k] = controls[k+N];
    forward_euler(agents[index].x, agents[index].y, agents[index].psi, agents[index].v_w, agents[index].psi_w, agents[index].u_phi, agents[index].u_v, Dt_MPC);
    
    //EVALUATE COST FUNCTION
    float F = 0;
    float total_cost = 0.0;
    
    float d[I];
    signed char B[Nb_x][Nb_y];
    //Accumulate the cells visited by the agents and the ones that the other agents plan to visit
    for (int i = 0; i < I; i++)
    {
        if (i != index) 
        {
            for (int j = 0; j < Nb_x; j++)
            {
                for (int k = 0; k < Nb_y; k++)
                {
                    B[j][k] = agents[index].B[j][k] || agents[i].B[j][k] || agents[i].b[j][k];
                }
            }
        }
    } 
    //anti colision for future steps. current step doesn't matter
    for (int k = 1; k < N; k++)
    {
      for (int i = 0; i < I; i++) 
      {
        d[i] = sqrt(pow(agents[index].x[k] - agents[i].x[k], 2) + pow(agents[index].y[k] - agents[i].y[k], 2));
        if ((d[i] < r_c) && (i != index)) return 99999999999999999999999999999999999999999999999999999999999999999999.9;
      }
    }
    for (int k = 0; k < (N+1); k++)
    { 
        //Lagrangian term
        
        //check if cells are being visited
        int x_idx = (int)agents[index].x[k]/res;
        int y_idx = (int)agents[index].y[k]/res;
        for (int i = (x_idx - R/res); i <= (x_idx + R/res); i++)
        {
            for (int j = (y_idx - R/res); j <= (y_idx + R/res); j++)
            {
                //printf("x=%f,y=%f \n",x1[k],y1[k]);
                if ( (i>=0) && (i<Nb_x) && (j>=0) && (j<Nb_y) && check_if_grid_is_inside_R(agents[index].x[k],agents[index].y[k],r_out[i][j]))
                {
                    if (B[i][j] == 0) B[i][j] = 1;
                    if (B[i][j] == 1) B[i][j] = -1;
                }
            }
        }
    }

    float phi_y;
    float sum_phi_y = 0;
    for (int k = 0; k < (N+1); k++)
    {
        phi_y = 0;
        for (int i = 0; i < Nb_x; i++)
        {
            for (int j = 0; j < Nb_y; j++)
            {
                if (B[i][j] != 0) phi_y = phi_y + phi_out[i][j]*B[i][j];
            }
        }
        sum_phi_y = sum_phi_y + phi_y;
    }
    F = 0.01*calculate_F(agents[index].x[N],agents[index].y[N],B); //r_x[Nb], r_y[Nb]: position of center of cells; b[Nb] is the binary of cell visit
    total_cost = F - sum_phi_y; //Eq 13
    return total_cost;
}
//ok
void forward_euler(float *x, float *y, float *psi, float v_w, float psi_w, float u_phi[], float u_v[], float Dt)
{   
    float dx, dy, v_g[N+1], chi[N+1], dchi;

    v_g[0] = sqrt( pow(u_v[0]*cos(psi[0])+v_w*cos(psi_w), 2) + pow(u_v[0]*sin(psi[0])+v_w*sin(psi_w), 2) );

    chi[0] = atan2( ( u_v[0]*sin(psi[0]) + v_w*sin(psi_w) ) , ( u_v[0]*cos(psi[0]) + v_w*cos(psi_w) ));
    if (chi[0] < 0) chi[0] = chi[0] + 2*pi;

  //FORWARD EULER METHOD TO INTEGRATE THE STATES OF AGENT (Eq. 02) (Eq. 04)
  for (int k = 0; k < N; k++)
  {
    //x
    dx = v_g[k]*cos(chi[k]);
    x[k+1] = x[k] + dx*Dt;
    //y
    dy = v_g[k]*sin(chi[k]);
    y[k+1] = y[k] + dy*Dt;
    //chi
    dchi = g*tan(u_phi[k])*cos(chi[k] - psi[k])/v_g[k];
    chi[k+1] = chi[k] + dchi*Dt;
    if (chi[k+1]>2*pi) chi[k+1] = chi[k+1]-2*pi;
    //psi
    psi[k+1] = chi[k+1] - asin(v_w/u_v[k+1]*sin(psi_w - chi[k+1]));
    //v_g
    v_g[k+1] = sqrt( pow(u_v[k+1]*cos(psi[k+1])+v_w*cos(psi_w), 2) + pow(u_v[k+1]*sin(psi[k+1])+v_w*sin(psi_w), 2) );
  }
}

// Get random between low and high
float getRandom(float low, float high)
{
    return (rand() % int(high*1000000 - low*1000000 + 1) + low*1000000)/1000000;
}

// Get random between 0.0f and 1.0f inclusive
float getRandomClamped()
{
    return (rand() % int(1*1000000-0*1000000+1*1000000) + 0*1000000)/1000000;
}

//Optimization
void PSO(int index, struct agent *agents)
{
    //update the cells that the agents are planning to visit
    for (int i = 0; i < I; i++)
    {
      if (i != index)
      {
        for (int k = 0; k < (N+1); k++)
        { 
          //check if cells are going to be visited
          int x_idx = (int)agents[i].x[k]/res;
          int y_idx = (int)agents[i].y[k]/res;
          for (int l = (x_idx - R/res); l <= (x_idx + R/res); l++)
          {
            for (int j = (y_idx - R/res); j <= (y_idx + R/res); j++)
            {
              //printf("x=%f,y=%f \n",x1[k],y1[k]);
              if ( (l>=0) && (l<Nb_x) && (j>=0) && (j<Nb_y) && check_if_grid_is_inside_R(agents[i].x[k],agents[i].y[k],r_out[l][j]))
              {
                agents[i].b[l][j] = 1;
              }
            }
          }
        }
      }
    }
      

    //restart the random
    srand(time(NULL));

    // Initialize particles
    for (int i = 0; i < NUM_OF_PARTICLES * NUM_OF_DIMENSIONS; i++)
    {
        //INITIALIZE THE POSITION OF THE PARTICLES
        if ((i % (N*2)) < N) positions[i] = getRandom(phi_min, phi_max);
        if ((i % (N*2)) >= N) positions[i] = getRandom(v_min, v_max);

        pBests[i] = positions[i];
        velocities[i] = 0;
    }
    //first particle is first gbest
    for (int i = 0; i < N*2; i++) gBest[i] = positions[i];

    for (int i = 0; i < NUM_OF_PARTICLES * NUM_OF_DIMENSIONS; i += NUM_OF_DIMENSIONS)
    {
        for(int k = 0; k < NUM_OF_DIMENSIONS; k++)
            temp[k] = pBests[i + k];
    
        if (cost_function(index, temp, agents) < cost_function(index, gBest, agents))
        {
            for (int k = 0; k < NUM_OF_DIMENSIONS; k++)
                gBest[k] = temp[k];
        }
    }

    cuda_pso(index, positions, velocities, pBests, gBest, agents, r_out, phi_out);

    for (int k = 0; k < N; k++) agents[index].u_phi[k] = gBest[k];
    for (int k = 0; k < N; k++) agents[index].u_v[k] = gBest[k+N];

    //zero all the cells that the agents are planning to visit
    for (int i = 0; i < I; i++)
    {
      for (int l = 0; l < Nb_x; l++)
      {
        for (int j = 0; j < Nb_y; j++)
        {
          agents[i].b[l][j] = 0;
        }
      }
    }
}

float calculate_F(float x,float y, signed char b_[Nb_x][Nb_y])
{
  float d = 0;
  for (int i = 0; i < Nb_x; i++)
  {
    for (int j = 0; j < Nb_y; j++)
    {
      if (b_[i][j] == 0) d = d + phi_out[i][j]*((x-r_out[i][j].x)*(x-r_out[i][j].x) + (y-r_out[i][j].y)*(y-r_out[i][j].y));
    }
  } 
  return d;
}

bool check_if_grid_is_inside_R(float x,float y, struct point2d r)
{
    float rx1 = r.x + res/2;
    float rx2 = r.x + res/2;
    float rx3 = r.x - res/2;
    float rx4 = r.x - res/2;

    float ry1 = r.y - res/2;
    float ry2 = r.y + res/2;
    float ry3 = r.y + res/2;
    float ry4 = r.y - res/2;

    if ( (calc_d(x,y,rx1,ry1) < (R2)) && (calc_d(x,y,rx2,ry2) < (R2)) && (calc_d(x,y,rx3,ry3) < (R2)) && (calc_d(x,y,rx4,ry4) < (R2)) )
    {
        return 1;
    }
    else 
    {
        return 0;
    }
}

float calc_d(float x1,float y1,float x2,float y2)
{
    return ((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
}