/*
 *  ou_generator.cpp
 *
 *  This file is part of NEST.
 *
 *  Copyright (C) 2004 The NEST Initiative
 *
 *  NEST is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  NEST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with NEST.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "ou_generator.h"

// Includes from libnestutil:
#include "logging.h"
#include "numerics.h"

// Includes from nestkernel:
#include "event_delivery_manager_impl.h"
#include "kernel_manager.h"
#include "universal_data_logger_impl.h"

// Includes from sli:
#include "dict.h"
#include "dictutils.h"
#include "doubledatum.h"
#include "integerdatum.h"

namespace nest
{
RecordablesMap< ou_generator > ou_generator::recordablesMap_;

template <>
void
RecordablesMap< ou_generator >::create()
{
  insert_( Name( names::I ), &ou_generator::get_I_avg_ );
}
}

/* ----------------------------------------------------------------
 * Default constructors defining default parameter
 * ---------------------------------------------------------------- */

nest::ou_generator::Parameters_::Parameters_()
  : mean_( 0.0 )    // pA
  , std_( 0.0 )     // pA / sqrt(s)
  // comment some variables, 重要 
  /*
  , std_mod_( 0.0 ) // pA / sqrt(s)
  , freq_( 0.0 )    // Hz
  , phi_deg_( 0.0 ) // degree
  */
  , dt_( Time::ms( 1.0 ) )
  , num_targets_( 0 )
  //added for tau, 重要
  , tau_( 0.0)      // tau
{
}

nest::ou_generator::Parameters_::Parameters_( const Parameters_& p )
  : mean_( p.mean_ )
  , std_( p.std_ )
  //, std_mod_( p.std_mod_ )
  , freq_( p.freq_ )
  , phi_deg_( p.phi_deg_ )
  , dt_( p.dt_ )
  , num_targets_( 0 ) 
  //added for tau, 重要
  , tau_( p.tau_ )// we do not copy connections
{
  // do not check validity of dt_ here, otherwise we cannot copy
  // to temporary in set(); see node copy c'tor
  dt_.calibrate();
}

nest::ou_generator::Parameters_& nest::ou_generator::Parameters_::
operator=( const Parameters_& p )
{
  if ( this == &p )
  {
    return *this;
  }

  mean_ = p.mean_;
  std_ = p.std_;
  /*
  std_mod_ = p.std_mod_;
  freq_ = p.freq_;
  phi_deg_ = p.phi_deg_;
  */
  dt_ = p.dt_;
  //added for tau, 重要
  tau_ = p.tau_;
  return *this;
}

nest::ou_generator::State_::State_()
  : y_0_( 0.0 )
  , y_1_( 0.0 )   // pA
  , I_avg_( 0.0 ) // pA
{
}

nest::ou_generator::Buffers_::Buffers_( ou_generator& n )
  : logger_( n )
{
}

nest::ou_generator::Buffers_::Buffers_( const Buffers_&, ou_generator& n )
  : logger_( n )
{
}

/* ----------------------------------------------------------------
 * Parameter extraction and manipulation functions
 * ---------------------------------------------------------------- */

void
nest::ou_generator::Parameters_::get( DictionaryDatum& d ) const
{
  ( *d )[ names::mean ] = mean_;
  ( *d )[ names::std ] = std_;
  // comment some variables, 重要   
  /*
  ( *d )[ names::std_mod ] = std_mod_;
  ( *d )[ names::dt ] = dt_.get_ms();
  ( *d )[ names::phase ] = phi_deg_;
  */
  ( *d )[ names::frequency ] = freq_;
  //added for tau, 重要
  ( *d )[ names::tau ] = tau_;
}

void
nest::ou_generator::State_::get( DictionaryDatum& d ) const
{
  ( *d )[ names::y_0 ] = y_0_;
  ( *d )[ names::y_1 ] = y_1_;
}

void
nest::ou_generator::Parameters_::set( const DictionaryDatum& d,
  const ou_generator& n )
{
  updateValue< double >( d, names::mean, mean_ );
  updateValue< double >( d, names::std, std_ );
  // comment some variables, 重要 
  /*
  updateValue< double >( d, names::std_mod, std_mod_ );
  updateValue< double >( d, names::frequency, freq_ );
  updateValue< double >( d, names::phase, phi_deg_ );
  */
  //added for tau, 重要
  updateValue< double >( d, names::tau, tau_ );
  double dt;
  if ( updateValue< double >( d, names::dt, dt ) )
  {
    dt_ = Time::ms( dt );
  }
  if ( std_ < 0 )
  {
    throw BadProperty( "The standard deviation cannot be negative." );
  }
  /*
  if ( std_mod_ < 0 )
  {
    throw BadProperty( "The standard deviation cannot be negative." );
  }
  
  if ( std_mod_ > std_ )
  {
    throw BadProperty(
      "The modulation apmlitude must be smaller or equal to the baseline "
      "amplitude." );
  }
  */
  if ( not dt_.is_step() )
  {
    throw StepMultipleRequired( n.get_name(), names::dt, dt_ );
  }
}


/* ----------------------------------------------------------------
 * Default and copy constructor for node
 * ---------------------------------------------------------------- */

nest::ou_generator::ou_generator()
  : Node()
  , device_()
  , P_()
  , S_()
  , B_( *this )
{
  recordablesMap_.create();
  if ( not P_.dt_.is_step() )
  {
    throw InvalidDefaultResolution( get_name(), names::dt, P_.dt_ );
  }
}

nest::ou_generator::ou_generator( const ou_generator& n )
  : Node( n )
  , device_( n.device_ )
  , P_( n.P_ )
  , S_( n.S_ )
  , B_( n.B_, *this )
{
  if ( not P_.dt_.is_step() )
  {
    throw InvalidTimeInModel( get_name(), names::dt, P_.dt_ );
  }
}


/* ----------------------------------------------------------------
 * Node initialization functions
 * ---------------------------------------------------------------- */

void
nest::ou_generator::init_state_( const Node& proto )
{
  const ou_generator& pr = downcast< ou_generator >( proto );

  device_.init_state( pr.device_ );
}

void
nest::ou_generator::init_buffers_()
{
  device_.init_buffers();
  B_.logger_.reset();

  B_.next_step_ = 0;
  B_.amps_.clear();
  B_.amps_.resize( P_.num_targets_, 0.0 );
}

void
nest::ou_generator::calibrate()
{
  B_.logger_.init();

  device_.calibrate();
  if ( P_.num_targets_ != B_.amps_.size() )
  {
    LOG( M_INFO,
      "ou_generator::calibrate()",
      "The number of targets has changed, drawing new amplitudes." );
    init_buffers_();
  }

  V_.dt_steps_ = P_.dt_.get_steps();

  const double h = Time::get_resolution().get_ms();
  const double t = kernel().simulation_manager.get_time().get_ms();
  
  // new added line for the e(-t/tau),重要
  S_.y_0_ = std::exp(-t/tau);
/* 
 * comment some variables, 重要
  // scale Hz to ms
  const double omega = 2.0 * numerics::pi * P_.freq_ / 1000.0;
  const double phi_rad = P_.phi_deg_ * 2.0 * numerics::pi / 360.0;

  // initial state
  S_.y_0_ = std::cos( omega * t + phi_rad );
  S_.y_1_ = std::sin( omega * t + phi_rad );

  // matrix elements
  V_.A_00_ = std::cos( omega * h );
  V_.A_01_ = -std::sin( omega * h );
  V_.A_10_ = std::sin( omega * h );
  V_.A_11_ = std::cos( omega * h );
   
*/  
}


/* ----------------------------------------------------------------
 * Update function and event hook
 * ---------------------------------------------------------------- */

nest::port
nest::ou_generator::send_test_event( Node& target,
  rport receptor_type,
  synindex syn_id,
  bool dummy_target )
{
  device_.enforce_single_syn_type( syn_id );

  if ( dummy_target )
  {
    DSCurrentEvent e;
    e.set_sender( *this );
    return target.handles_test_event( e, receptor_type );
  }
  else
  {
    CurrentEvent e;
    e.set_sender( *this );
    const port p = target.handles_test_event( e, receptor_type );
    if ( p != invalid_port_ and not is_model_prototype() )
    {
      ++P_.num_targets_;
    }
    return p;
  }
}

//
// Time Evolution Operator
//
void
nest::ou_generator::update( Time const& origin,
  const long from,
  const long to )
{
  assert(
    to >= 0 && ( delay ) from < kernel().connection_manager.get_min_delay() );
  assert( from < to );

  const long start = origin.get_steps();

  for ( long offs = from; offs < to; ++offs )
  {
    S_.I_avg_ = 0.0;

    const long now = start + offs;

    if ( not device_.is_active( Time::step( now ) ) )
    {
      B_.logger_.record_data( origin.get_steps() + offs );
      continue;
    }

    // >= in case we woke from inactivity
    if ( now >= B_.next_step_ )
    {
      // compute new currents
      for ( AmpVec_::iterator it = B_.amps_.begin(); it != B_.amps_.end();
            ++it )
      {
		// set the initial value,重要
        *(B_.amps_.begin()) = 0;
        // Changed the formula to compute the currents according to the equation 4, page 3,重要
        *it = *(it-1) * S_.y_
          + P_.std_ * std::sqrt( 1 - S_.y_0_* S_.y_0_ ) / std::sqrt(2*tau)
            * V_.normal_dev_( kernel().rng_manager.get_rng( get_thread() ) );
        S_.I_avg_ += *it;
      }
      S_.I_avg_ /= B_.amps_.size();
      B_.logger_.record_data( origin.get_steps() + offs );

      // use now as reference, in case we woke up from inactive period
      B_.next_step_ = now + V_.dt_steps_;
    }

    DSCurrentEvent ce;
    kernel().event_delivery_manager.send( *this, ce, offs );
  }
}

void
nest::ou_generator::event_hook( DSCurrentEvent& e )
{
  // get port number
  const port prt = e.get_port();

  // we handle only one port here, get reference to vector elem
  assert( 0 <= prt && static_cast< size_t >( prt ) < B_.amps_.size() );

  e.set_current( B_.amps_[ prt ] );
  e.get_receiver().handle( e );
}

void
nest::ou_generator::handle( DataLoggingRequest& e )
{
  B_.logger_.handle( e );
}
