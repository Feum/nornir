# ---------------------------------------------------------------------------
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License version 2 as 
#  published by the Free Software Foundation.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
#
#  As a special exception, you may use this file as part of a free software
#  library without restriction.  Specifically, if other files instantiate
#  templates or use macros or inline functions from this file, or you compile
#  this file and link it with other files to produce an executable, this
#  file does not by itself cause the resulting executable to be covered by
#  the GNU General Public License.  This exception does not however
#  invalidate any other reasons why the executable file might be covered by
#  the GNU General Public License.
#
# ---------------------------------------------------------------------------


CXXFLAGS += -Wall `pkg-config --cflags opencv` #-DTRACE_FASTFLOW  -pg
OPENCVLIB = `pkg-config --libs opencv` -pthread 
LDLIBS 	+= $(OPENCVLIB)

TARGETS	= seqvideo seqvideo_pedestrian video video_nocvthreads 

all: $(TARGETS)

seqvideo: seqvideo.cpp
	$(CXX) $(OPTIMIZE_FLAGS) $(CXXFLAGS) -o $@ $< $(INCS) $(LDFLAGS) $(LDLIBS)
seqvideo_pedestrian: seqvideo_pedestrian.cpp
	$(CXX) $(OPTIMIZE_FLAGS) $(CXXFLAGS) -DSHOWTIMES -DNO_CV_THREADS -o $@ $< $(INCS) $(LDFLAGS) $(LDLIBS)
video: video.cpp
	$(CXX) -DNO_DEFAULT_MAPPING $(OPTIMIZE_FLAGS) $(CXXFLAGS) -o $@ $< $(INCS) $(LDFLAGS) $(LDLIBS)
video_nocvthreads: video.cpp
	$(CXX) -DNO_DEFAULT_MAPPING -DNO_CV_THREADS $(OPTIMIZE_FLAGS) $(CXXFLAGS) -o $@ $< $(INCS) $(LDFLAGS) $(LDLIBS)
clean: 
	rm -f $(TARGETS)
cleanall: clean
	rm -f *.o *.a *~ 

