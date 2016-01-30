/* -*- Mode: C++; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: nil -*- */
/* ***************************************************************************
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License version 2 as 
 *  published by the Free Software Foundation.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 *
 *  As a special exception, you may use this file as part of a free software
 *  library without restriction.  Specifically, if other files instantiate
 *  templates or use macros or inline functions from this file, or you compile
 *  this file and link it with other files to produce an executable, this
 *  file does not by itself cause the resulting executable to be covered by
 *  the GNU General Public License.  This exception does not however
 *  invalidate any other reasons why the executable file might be covered by
 *  the GNU General Public License.
 *
 ****************************************************************************
 */

/* 
 * Author: Marco Danelutto <marcod@di.unipi.it> 
 * Date:   September 2015
 * 
 */
//  Version using only the ordered farm:
//    ofarm(Stage1+Stage2)

#include <opencv2/opencv.hpp>
#include "../../manager.hpp"

using namespace ff; 
using namespace cv;
using namespace adpff;

// reads frame and sends them to the next stage
struct Source : AdaptiveNode {
    const std::string filename;
    VideoCapture cap;
    
    Source(const std::string filename):filename(filename) {
        cap = VideoCapture(filename.c_str());
        if(!cap.isOpened())  {
            std::cout << "Error opening input file" << std::endl;
            exit(-1);
        }    
    }
  
    
    void * svc(void *) {
        Mat * frame = new Mat();
        if(!cap.read(*frame)){
            std::cout << "End of stream in input" << std::endl; 
            TERMINATE_APPLICATION;
        }else{
            return (void*) frame;
        }
    }
}; 

// this stage applys all the filters:  the GaussianBlur filter and the Sobel one, 
// and it then sends the result to the next stage
struct Stage1 : AdaptiveNode {
    void * svc(void *task) {
        Mat* frame = (Mat*) task;
        Mat frame1;
        //cv::GaussianBlur(*frame, frame1, cv::Size(0, 0), 3);
        //cv::addWeighted(*frame, 1.5, frame1, -0.5, 0, *frame);
        cv::Sobel(*frame,*frame,-1,1,0,3);
        return (void*) frame;
    }
}; 

// this stage shows the output
struct Drain: AdaptiveNode {
    Drain(bool ovf):outvideo(ovf) {}

    int svc_init() {
	if(outvideo) namedWindow("edges",1);
	return 0; 
    }

    void *svc (void * task) {
        Mat* frame = (Mat*) task;
	if(outvideo) {
	    imshow("edges", *frame);
	    waitKey(30);    
	} 
	delete frame;
	return (void*) GO_ON;
    }
protected:
    const bool outvideo; 
}; 

int main(int argc, char *argv[]) {
    //ffvideo input.mp4 filterno output nw1
    Mat edges;

    if(argc == 1) {
      std::cout << "Usage is: " << argv[0] 
                << " input_filename videooutput nw1" 
		<< std::endl; 
      return(0); 
    }
    
    // output 
    bool outvideo = false; 
    if(atoi(argv[2]) == 1) outvideo = true; 
    
    // pardegree 
    size_t nw1 = 1;
    if(argc == 4) {
      nw1 = atol(argv[3]); 
    }

    // creates an ordered farm
#if 0
    ff_OFarm<cv::Mat> ofarm( [nw1]() {
            
            std::vector<std::unique_ptr<ff_node> > W; 
            for(size_t i=0; i<nw1; i++) 
                W.push_back(make_unique<Stage1>());
            return W;
            
        } ());
#endif
    ff_ofarm ofarm;
    std::vector<ff_node*> W;
    for(size_t i=0; i<nw1; i++)
        W.push_back(new Stage1());
    ofarm.add_workers(W);
    
    Source source(argv[1]);
    ofarm.setEmitterF((ff_node*) &source);
    Drain  drain(outvideo);
    ofarm.setCollectorF((ff_node*) &drain);
    
    adpff::Observer obs;
    adpff::Parameters ap("parameters.xml", "archdata.xml");
    ap.observer = &obs;
    adpff::ManagerFarm<ofarm_lb, ofarm_gt> amf(&ofarm, ap);

    amf.start();
    amf.join();

    return 0;
}

