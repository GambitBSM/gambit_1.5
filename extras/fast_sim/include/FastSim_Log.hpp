#ifndef LOG_HPP
#define LOG_HPP

// very simple logger that is thread safe provided as an example in an issue of Dr Dobbs
// Class and template denitifions - does require C++11
// A. Saavedra  June 25-06-2014
#include <iostream>
#include <string>
#include <sstream>
#include <memory>
#include <mutex>

namespace logging
{

	enum severity_type
	{
    info = 1,
		debug,
		error,
		warning
	};


  /*
	class log_policy_interface
	{
	public:
		virtual void		open_ostream(const std::string& name) = 0;
		virtual void		close_ostream() = 0;
		virtual void		write(const std::string& msg) = 0;

	};
  

	// *
	// * Implementation which allow to write into a file
	// *

	class file_log_policy : public log_policy_interface
	{
		std::unique_ptr< std::ofstream > out_stream;
	public:
		file_log_policy() : out_stream( new std::ofstream ) {}
		void open_ostream(const std::string& name);
		void close_ostream();
		void write(const std::string& msg);
		~file_log_policy();
	};
  */


	//template< typename log_policy >
	class logger
	{
		unsigned _log_line_number;
		std::string get_time();
		std::string get_logline_header();
		//std::stringstream _log_stream;
    int _verbosity_level;
		std::mutex _write_mutex;

		//Core printing functionality
		void print_impl();
		template<typename First, typename...Rest>
		void print_impl(First parm1, Rest...parm);
	public:
		logger(int verbosity );
		logger();

    void set_verbosity(int verbosity);
		template< severity_type severity,int verbosity, typename...Args >
		void print(Args...args );

		~logger();
	};


  
		template< severity_type severity, int verbosity, typename...Args >
	void logger::print(Args...args )
	{
    std::string messg;
		_write_mutex.lock();

    switch( severity )
    {
      case severity_type::debug:
        messg = "<FastSim DEBUG> :";
        break;
      case severity_type::info:
        messg = "<FastSim INFO> :";
        break;
      case severity_type::warning:
        messg = "<FastSim WARNING> :";
        break;
      case severity_type::error:
        messg ="<FastSim ERROR> :";
        break;
    }

    if (severity ==severity_type::debug) {
      if (verbosity <= _verbosity_level) {
        print_impl( messg,args... );
      }
    }
    else {
      print_impl( messg,args... );
    }
    
    _write_mutex.unlock();
  }

		template<typename First, typename...Rest >
	void logger::print_impl(First parm1, Rest...parm)
	{
    std::cout <<" "<<parm1;
		print_impl(parm...);	
	}



}
#endif

