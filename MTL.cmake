#
# Copyright (c) 2021 Seung Jae Lee
#

# Change the default for message().
function(message)
  if(${ARGC} EQUAL 1)
    _message(STATUS ${ARGV})
  else()
    _message(${ARGV})
  endif()
endfunction()

message("Applied MTL.cmake")
