

# https://midas3.kitware.com/midas/download/item/318227/castxml-linux.tar.gz
# https://midas3.kitware.com/midas/download/item/318762/castxml-macosx.tar.gz


# Arguments:
#   1. BOSS directory
#   2. cmake command
#   3. primary URL


if [ -d $1/castxml ] ; then
  echo "DEBUG: Found the castxml directory"
else
  echo "DEBUG: Did NOT find the castxml directory. Will try to download castxml now"

	# Download
	echo "DEBUG: arg 1: $1"
	echo "DEBUG: arg 2: $2"
	echo "DEBUG: arg 3: $3"
	axel_worked=0
	filename=$($2 -E echo $3 | sed 's#.*/##g')
	$2 -E make_directory $1 >/dev/null
	# Go to wget/curl if axel is not present
	if command -v axel >/dev/null; then
	  # Go to wget/curl if POST data have been provided
	  if [ -z "$9" ]; then
	    if $2 -E chdir $1 axel $3; then
	      axel_worked=1
	    else
	      $2 -E echo "Axel failed! The link probably redirects to https. Falling back to wget/curl..."
	    fi
	  fi
	fi
	if [ "${axel_worked}" = "0" ]; then
	  if command -v wget >/dev/null; then
	    if [ -z "$9" ]; then
	      wget $3 -O $1/${filename}
	    else
	      wget --post-data "$9" ${10} -O $1/${filename}
	    fi
	  elif command -v curl >/dev/null; then
	    if [ -z "$9" ]; then
	      $2 -E chdir $1 curl -O $3
	    else
	      $2 -E chdir $1 curl -O -c $cfile --data "$9" ${10}
	      $2 -E chdir $1 curl -O -b $cfile $3
	      $2 -E remove $1/$cfile
	    fi
	  else
	    $2 -E cmake_echo_color --red --bold "ERROR: No axel, no wget, no curl?  What kind of OS are you running anyway?"
	    exit 1
	  fi
	fi

fi

# Do the extraction
$2 -E echo "Extracting $1/${filename}"
cd $1
$2 -E tar -xf $1/${filename}

