#include <opencv2/opencv.hpp>
#include <stdint.h>
#include "input_videos/video_data.h"

using namespace cv;

int main()
{
    Mat output;

    // erstes Frame
    Mat prev(VIDEO_H, VIDEO_W, CV_8UC1, (void*)video_data[0]);
    cvtColor(prev, output, COLOR_GRAY2BGR);

    for(int t = 1; t < VIDEO_T; t++)
    {
        Mat next(VIDEO_H, VIDEO_W, CV_8UC1, (void*)video_data[t]);

        std::vector<Point2f> pointsPrev;
        goodFeaturesToTrack(prev, pointsPrev, 300, 0.01, 10);

        std::vector<Point2f> pointsNext;
        std::vector<uchar> status;
        std::vector<float> err;

        calcOpticalFlowPyrLK(prev, next, pointsPrev, pointsNext, status, err);

        for(size_t i = 0; i < pointsPrev.size(); i++)
        {
            if(status[i])
            {
                line(output, pointsPrev[i], pointsNext[i], Scalar(0,255,0), 2);
                circle(output, pointsNext[i], 2, Scalar(0,0,255), -1);
            }
        }

        prev = next;
    }

    imwrite("optical_flow.png", output);

    return 0;
}
