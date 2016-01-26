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

#include <ff/utils.hpp>
#include <opencv2/opencv.hpp>
#include <mammut/utils.hpp>


// #define SHOWTIMES
#ifdef SHOWTIMES
/* times are in millisecond */
#define TIME(t0)  ( (ff::getusec() - t0) /  1000 )
#endif

using namespace cv;

int main(int argc, char *argv[]) {
    
    if(argc == 1) {
        std::cout << "Usage is: " << argv[0]
              << " input_filename filterno videooutput"
              << std::endl;
        return(0);
    }
        
    // input file
    //VideoCapture cap(1); // open the default camera
    VideoCapture cap(argv[1]);

    cv::VideoWriter outputFile;

    if(!cap.isOpened())  {  // check if we succeeded
        std::cerr << "Error opening input file" << std::endl;
        return -1;
    }
    std::cout << "Input file " << argv[1] << " opened" << std::endl;
    
    // filter parameters
    bool filter1=false, filter2=false;
    if(atoi(argv[2]) == 1) { 
        filter1 = true;  filter2 = false; 
        std::cout << "Applying enhnace filter only" << std::endl;
    }
    if(atoi(argv[2]) == 2) { 
        filter1 = false; filter2 = true;     
        std::cout << "Applying emboss filter only" << std::endl;
    }
    if(atoi(argv[2]) == 3) { 
        filter1 = true;  filter2 = true;     
        std::cout << "Applying both filters" << std::endl;
    }
    
    
    Mat edges;
    // output 
    switch(atoi(argv[3])){
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

    cv::HOGDescriptor hog;
    hog.setSVMDetector(cv::HOGDescriptor::getDefaultPeopleDetector());

    
    ff::ffTime(ff::START_TIME);
    int frames = 0; 
    for(;;)  {
        Mat frame1;
        Mat frame;
        
#ifdef SHOWTIMES
        unsigned long t0 = ff::getusec();
#endif
        if(cap.read(frame) == false) 
            break; 
#ifdef SHOWTIMES
        std::cout << "Read " << TIME(t0) << std::endl;
#endif
        
        frames++; 

        // start and end times
        double start, end;

        // fps calculated using number of frames / seconds
        double fps;

        // frame counter
        int counter = 0;

        // floating point seconds elapsed since start
        double sec;

        // start the clock
        start = mammut::utils::getMillisecondsTime();
        vector<Rect> found, found_filtered;
        hog.detectMultiScale(frame, found, 0, Size(8,8), Size(32,32), 1.05, 2);
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
            rectangle(frame, r.tl(), r.br(), Scalar(0,255,0), 3);
        }

        // see how much time has elapsed
        end = mammut::utils::getMillisecondsTime();

        // calculate current FPS
        ++counter;
        sec = (end - start) / 1000.0;

        fps = counter / sec;

        // will print out Inf until sec is greater than 0
        printf("FPS = %.2f Sec = %.2f\n", fps, sec);

#ifdef SHOWTIMES
        t0 = ff::getusec();
#endif
        switch(atoi(argv[3])){
        case 0:{
            ;
        }break;
        case 1:{
            imshow("edges", frame);
            if(waitKey(30) >= 0) break;
        }break;
        case 2:{
            outputFile.write(frame);
        }break;
        }
#ifdef SHOWTIMES
        std::cout << "Show " << TIME(t0) << std::endl;
#endif
    }
    ff::ffTime(ff::STOP_TIME);
    std::cout << "Elapsed time is " << ff::ffTime(ff::GET_TIME) << " ms\n";
    std::cout << "Average time per frame " << (ff::ffTime(ff::GET_TIME) / frames) <<  " ms\n";
    std::cout << "(with " << frames << " frames)" << std::endl;
    
    // the camera will be deinitialized automatically in VideoCapture destructor
    return 0;
}
