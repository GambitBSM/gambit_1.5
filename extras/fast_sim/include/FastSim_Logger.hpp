#ifndef LOGGER_HPP
#define LOGGER_HPP

// very simple logger that is thread safe provided as an example in an issue of Dr Dobbs
// Macro definitions and instantiation of the logger - does require C++11
// A. Saavedra  June 25-06-2014


#include <iostream>
#include "FastSim_Log.hpp"
#define LOGGING_LEVEL_1 1
//static logging::logger log_inst( "fastsim_execution.log" );

#ifdef LOGGING_LEVEL_1
#define LOG_INFO log_inst.print< logging::severity_type::info,0 >
#define LOG_DEBUG1 log_inst.print<logging::severity_type::debug,1 >
//#define LOG_DEBUG1 log_inst.print<logging::severity_type::debug >
#define LOG_DEBUG2 log_inst.print<logging::severity_type::debug,2 >
#define LOG_ERR log_inst.print< logging::severity_type::error,0 >
#define LOG_WARN log_inst.print< logging::severity_type::warning,0 >
#else
#define LOG_INFO log_inst.print< logging::severity_type::info,0 >
#define LOG_DEBUG(...) 
#define LOG_ERR(...)
#define LOG_WARN(...)
#endif

/*
#ifdef LOGGING_LEVEL_2
#define ELOG log_inst.print< logging::severity_type::debug >
#define ELOG_ERR log_inst.print< logging::severity_type::error >
#define ELOG_WARN log_inst.print< logging::severity_type::warning >
#else
#define ELOG(...) 
#define ELOG_ERR(...)
#define ELOG_WARN(...)
#endif
*/

#endif
