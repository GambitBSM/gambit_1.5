# - Finds RestFrames instalation
# This module sets up RestFrames information
# It defines:
# ...
set(ver "1.0.1")
set(dir "${PROJECT_SOURCE_DIR}/contrib/RestFrames-${ver}")
unset(RestFrames_1.0.1_LIBRARY CACHE)
find_library(RestFrames_${ver}_LIBRARY RestFrames 
  ${dir}/lib/
)
