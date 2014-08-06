#include "FastSim_Log.hpp"

namespace logging
{

  /*
  void file_log_policy::open_ostream(const std::string& name)
	{
		out_stream->open( name.c_str(), std::ios_base::binary|std::ios_base::out );
		if( !out_stream->is_open() ) 
		{
			throw(std::runtime_error("LOGGER: Unable to open an output stream"));
		}
	}

	void file_log_policy::close_ostream()
	{
		if( out_stream )
		{
			out_stream->close();
		}
	}

	void file_log_policy::write(const std::string& msg)
	{
		(*out_stream)<<msg<<std::endl;
    // need to fix this fudge at some stage
    std::cout<<msg<<std::endl;
	}

	file_log_policy::~file_log_policy()
	{
		if( out_stream )
		{
			close_ostream();
		}
	}
  */

	void logger::print_impl()
	{
		//policy->write( get_logline_header() + log_stream.str() );
    //std::cout << get_logline_header() + _log_stream.str() << std::endl;
		//_log_stream.str("");
	}


	//template< typename log_policy >
	std::string logger::get_time()
	{
		std::string time_str;
		time_t raw_time;
		
		time( & raw_time );
		time_str = ctime( &raw_time );

		//without the newline character
		return time_str.substr( 0 , time_str.size() - 1 );
	}


//	template< typename log_policy >
	std::string logger::get_logline_header()
	{
		std::stringstream header;

		header.str("");
		header.fill('0');
		header.width(4);
		header << _log_line_number++ <<" < "<<get_time();

		//header.fill('0');
		//header.width(7);
		//header <<clock()<<" > ~ ";
		header <<" > ";

		return header.str();
	}


	//template< typename log_policy >
	logger::logger()
	{
		_log_line_number = 0;
    _verbosity_level = 0; // set by the user

//		policy = new log_policy;
//		if( !policy )
//		{
//			throw std::runtime_error("LOGGER: Unable to create the logger instance"); 
//		}

//		log_stream.precision(3);
//		policy->open_ostream( name );

	}



	logger::logger(int verbosity )
	{
		_log_line_number = 0;
    _verbosity_level = verbosity; // set by the user

//		policy = new log_policy;
//		if( !policy )
//		{
//			throw std::runtime_error("LOGGER: Unable to create the logger instance"); 
//		}

//		log_stream.precision(3);
//		policy->open_ostream( name );

	}

  void logger::set_verbosity(int verbosity)
  {
    _verbosity_level = verbosity;
  }

	//template< typename log_policy >
	logger::~logger()
	{

//		if( policy )
//		{
//			policy->close_ostream();
//			delete policy;
//		}
	}




}
