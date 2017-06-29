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
using namespace nornir;

// reads frame and sends them to the next stage
struct Source : AdaptiveNode {
    std::vector<std::string> filenames;
    VideoCapture* cap;
    int maxFrames;
    size_t currentFileId;
    size_t currentFrameId;


    void openFile(uint id){
        cap = new VideoCapture(filenames.at(id).c_str());
        if(!cap->isOpened())  {
            std::cout << "Error opening input file" << std::endl;
            exit(-1);
        }
    }

    Source(std::vector<std::string> filenames, int maxFramesPerFile):filenames(filenames), maxFrames(maxFramesPerFile){
        openFile(0);
        currentFileId = 0;
        currentFrameId = 0;
    }
  
    
    void * svc(void *) {
        Mat * frame = new Mat();
        if(!cap->read(*frame) || (maxFrames && currentFrameId >= (uint) maxFrames)){
            std::cout << "End of stream in input" << std::endl; 
            if(currentFileId + 1 < filenames.size()){
                ++currentFileId;
                cap->release();
                delete cap;
                openFile(currentFileId);
                currentFrameId = 0;
                cap->read(*frame);
            }else{
                TERMINATE_APPLICATION;
            }
        }

        ++currentFrameId;
        return (void*) frame;

    }
}; 

// this stage applys all the filters:  the GaussianBlur filter and the Sobel one, 
// and it then sends the result to the next stage
struct Stage1 : AdaptiveNode {
    double d;
    uint id;
    Stage1(double d, uint id):d(d),id(id){;}


    void * svc(void *task) {
        Mat* frame = (Mat*) task;
        Mat frame1;
        //cv::GaussianBlur(*frame, frame1, cv::Size(0, 0), 3);
        //cv::addWeighted(*frame, 1.5, frame1, -0.5, 0, *frame);
        //cv::Sobel(*frame,*frame,-1,1,0,3);

        cv::bilateralFilter(*frame , frame1, d, 80, 80);
        *frame = frame1;

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
    //ffvideo numframes d output nw1 file1 file2 ... filen

#ifdef NO_CV_THREADS
    setNumThreads(0);
#endif

    if(argc == 1) {
      std::cout << "Usage is: " << argv[0] 
                << " numframes(0 is all frames) d output nw1 file1 file2 ... filen"
		<< std::endl; 
      return(0); 
    }
    
    // output 
    bool outvideo = false; 
    int numframes = atoi(argv[1]);
    int d = atoi(argv[2]);
    if(atoi(argv[3]) == 1) 
        outvideo = true; 
    
    // pardegree 
    size_t nw1 = atol(argv[4]); 

    // creates an ordered farm
    ff_ofarm ofarm;
    std::vector<ff_node*> W;
    for(size_t i=0; i<nw1; i++)
        W.push_back(new Stage1(d, i+1));
    ofarm.add_workers(W);

    std::vector<std::string> filenames;
    uint numfiles = argc - 5;
    for(uint i = 0; i < numfiles; i++){
        filenames.push_back(argv[i + 5]);
    }

    Source source(filenames, numframes);
    ofarm.setEmitterF((ff_node*) &source);
    Drain  drain(outvideo);
    ofarm.setCollectorF((ff_node*) &drain);
    
    nornir::Parameters ap("parameters.xml");
    if(numframes){
        ap.expectedTasksNumber = numframes*numfiles;
    }
    nornir::ManagerFarm<ofarm_lb, ofarm_gt> amf(&ofarm, ap);

    amf.start();
    amf.join();

    return 0;
}

