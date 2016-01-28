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

#include <opencv2/opencv.hpp>
#include "../../src/manager.hpp"

using namespace ff; 
using namespace cv;
using namespace std;

// reads frame and sends them to the next stage
struct Source : adpff::AdaptiveNode {
    VideoCapture& _cap;
    uint _numTasks;
    Source(VideoCapture& cap):_cap(cap), _numTasks(0){}
  
    void* svc(void*) {
        Mat * frame = new Mat();
        if(_cap.read(*frame)){
            ++_numTasks;
            return (void*) frame;
        }

        std::cout << "End of stream in input" << std::endl;
        std::cout << _numTasks << " frames sent." << std::endl;
        TERMINATE_APPLICATION;
    }
}; 

// this stage applys all the filters:  the GaussianBlur filter and the Sobel one, 
// and it then sends the result to the next stage
struct Stage1 : adpff::AdaptiveNode {
    Stage1():nframe(0){
        hog.setSVMDetector(cv::HOGDescriptor::getDefaultPeopleDetector());
    }

    void * svc(void* task) {
        cv::Mat *frame = (cv::Mat*) task;
        vector<Rect> found, found_filtered;

        hog.detectMultiScale(*frame, found, 0, Size(8,8), Size(32,32), 1.05, 2);
        size_t i, j;
        for (i=0; i<found.size(); i++)
        {
            Rect r = found[i];
            for (j=0; j<found.size(); j++)
                if (j!=i && (r & found[j]) == r)
                    break;
            if (j== found.size())
                found_filtered.push_back(r);
        }

        for (i=0; i<found_filtered.size(); i++)
        {
            Rect r = found_filtered[i];
            r.x += cvRound(r.width*0.1);
            r.width = cvRound(r.width*0.8);
            r.y += cvRound(r.height*0.07);
            r.height = cvRound(r.height*0.8);
            rectangle(*frame, r.tl(), r.br(), Scalar(0,255,0), 3);
        }
        return (void*) frame;
    }
    long nframe;
    cv::HOGDescriptor hog;

}; 

// this stage shows the output
struct Drain: adpff::AdaptiveNode {
    Drain(uint out, cv::VideoWriter& outputFile):_out(out), _outputFile(outputFile){
        ;
    }

    void *svc (void* task) {
        cv::Mat * frame = (cv::Mat*) task;
        switch(_out){
        case 0:{
            ;
        }break;
        case 1:{
            cv::imshow("edges", *frame);
            if(waitKey(30) >= 0) break;
        }break;
        case 2:{
            _outputFile.write(*frame);
        }break;
        }
        delete frame;
        return GO_ON;
    }
protected:
    const uint _out;
    cv::VideoWriter& _outputFile;
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
    

    VideoCapture cap(argv[1]);
    if(!cap.isOpened())  {
        std::cout << "Error opening input file" << std::endl;
        return -1;
    }

    cv::VideoWriter outputFile;
    // output 
    switch(atoi(argv[2])){
    case 0:{
        ;
    }break;
    case 1:{
        namedWindow("edges",1);
    }break;
    case 2:{
        outputFile.open("output.mp4", cap.get(CV_CAP_PROP_FOURCC), cap.get(CV_CAP_PROP_FPS), cvSize((int)cap.get(CV_CAP_PROP_FRAME_WIDTH),(int)cap.get(CV_CAP_PROP_FRAME_HEIGHT)), true);
    }break;
    }
    
    // pardegree 
    size_t nw1 = 1;
    if(argc == 4) {
      nw1 = atol(argv[3]); 
    }

    // creates an ordered farm
    ff_OFarm<cv::Mat> ofarm( [nw1]() {
            
            std::vector<std::unique_ptr<ff_node> > W; 
            for(size_t i=0; i<nw1; i++) 
                W.push_back(make_unique<Stage1>());
            return W;
            
        } ());

    Source source(cap);
    ofarm.setEmitterF(source);
    Drain  drain(atoi(argv[2]), outputFile);
    ofarm.setCollectorF(drain);

    adpff::Observer obs;
    adpff::Parameters ap("parameters.xml", "archdata.xml");
    ap.observer = &obs;
    adpff::ManagerFarm<ofarm_lb, ofarm_gt> amf(&ofarm, ap);

    ffTime(START_TIME);
    amf.start();
    amf.join();
#if 0
    if (ofarm.run_and_wait_end()<0) {
        error("running farm");
        return -1;
    }
#endif
    ffTime(STOP_TIME);

    std::cout << "Elapsed (farm(" << nw1 << "): elapsed time =" ;     
    std::cout << ffTime(GET_TIME) << " ms\n";    
    return 0;
}

