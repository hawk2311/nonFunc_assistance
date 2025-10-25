#include <opencv4/opencv2/opencv.hpp>
#include <iostream>

//Source: https://chatgpt.com/share/68fd0b98-3244-8005-92da-f5dad4a98394
using namespace cv;
using namespace std;

// Simple program that reads a video, computes optical flow (Farneback by default),
// visualizes it, and writes an output video. Designed to be lightweight and portable.
// Build with: `g++ optical_flow_riscv.cpp `pkg-config --cflags --libs opencv4` -O2 -o optical_flow`

static void drawOpticalFlowArrows(const Mat& flow, Mat& out, int step=16, Scalar color=Scalar(0,255,0)){
    // Draw arrows on out image at each 'step' pixels using flow vectors
    for(int y = 0; y < flow.rows; y += step){
        for(int x = 0; x < flow.cols; x += step){
            const Point2f& f = flow.at<Point2f>(y, x);
            Point p(x, y);
            Point q(cvRound(x + f.x), cvRound(y + f.y));
            line(out, p, q, color, 1, LINE_AA);
            // draw arrow head
            double angle = atan2((double)p.y - q.y, (double)p.x - q.x);
            const double PI = 3.14159265358979323846;
            Point p1 = Point(cvRound(q.x + 4 * cos(angle + PI/6)), cvRound(q.y + 4 * sin(angle + PI/6)));
            Point p2 = Point(cvRound(q.x + 4 * cos(angle - PI/6)), cvRound(q.y + 4 * sin(angle - PI/6)));
            line(out, q, p1, color, 1, LINE_AA);
            line(out, q, p2, color, 1, LINE_AA);
        }
    }
}

static Mat flowToColor(const Mat& flow){
    // Convert 2-channel flow to HSV color visualization
    Mat hsv(flow.size(), CV_8UC3);
    for(int y=0;y<flow.rows;++y){
        for(int x=0;x<flow.cols;++x){
            Point2f f = flow.at<Point2f>(y,x);
            float fx = f.x;
            float fy = f.y;
            float magnitude = sqrt(fx*fx + fy*fy);
            float angle = atan2(fy, fx); // -pi..pi
            float hue = (angle + (float)CV_PI) * (180.0f / (2.0f*(float)CV_PI)); // 0..180
            float sat = 255;
            float val = std::min(255.0f, magnitude * 10.0f); // scale for visibility
            hsv.at<Vec3b>(y,x) = Vec3b((uchar)hue, (uchar)sat, (uchar)val);
        }
    }
    Mat bgr;
    cvtColor(hsv, bgr, COLOR_HSV2BGR);
    return bgr;
}

int main(int argc, char** argv){
    if(argc < 3){
        cerr << "Usage: " << argv[0] << " <input-video> <output-video> [--method farneback|lk]" << endl;
        return 1;
    }

    string inPath = argv[1];
    string outPath = argv[2];
    string method = "farneback";
    if(argc >= 4) method = argv[3];

    VideoCapture cap(inPath);
    if(!cap.isOpened()){
        cerr << "ERROR: could not open input video: " << inPath << endl;
        return 1;
    }

    int width = static_cast<int>(cap.get(CAP_PROP_FRAME_WIDTH));
    int height = static_cast<int>(cap.get(CAP_PROP_FRAME_HEIGHT));
    double fps = cap.get(CAP_PROP_FPS);
    if(fps <= 0) fps = 25; // fallback

    VideoWriter writer;
    int fourcc = VideoWriter::fourcc('a','v','c','1'); // try h264; on some platforms use 'M','J','P','G'
    bool opened = writer.open(outPath, fourcc, fps, Size(width, height));
    if(!opened){
        // fallback to MJPG
        fourcc = VideoWriter::fourcc('M','J','P','G');
        opened = writer.open(outPath, fourcc, fps, Size(width, height));
        if(!opened){
            cerr << "ERROR: could not open output video: " << outPath << endl;
            return 1;
        }
    }

    Mat prevGray, gray, frame;
    if(!cap.read(frame)){
        cerr << "ERROR: empty video or cannot read first frame" << endl;
        return 1;
    }
    cvtColor(frame, prevGray, COLOR_BGR2GRAY);

    // For LK tracking (sparse), we maintain points
    vector<Point2f> prevPts, nextPts;
    vector<unsigned char> status;
    vector<float> err;

    if(method == "lk"){
        // detect good features to track in the first frame
        goodFeaturesToTrack(prevGray, prevPts, 1000, 0.01, 7);
    }

    int frameIdx = 1;
    while(true){
        if(!cap.read(frame)) break;
        cvtColor(frame, gray, COLOR_BGR2GRAY);

        Mat display = frame.clone();

        if(method == "farneback"){
            Mat flow;
            calcOpticalFlowFarneback(prevGray, gray, flow,
                                     0.5, // pyr_scale
                                     3,   // levels
                                     15,  // winsize
                                     3,   // iterations
                                     5,   // poly_n
                                     1.2, // poly_sigma
                                     0);
            // dense visualization
            Mat flowColor = flowToColor(flow);
            addWeighted(display, 0.6, flowColor, 0.4, 0, display);
            // optionally draw arrows every 16 pixels
            drawOpticalFlowArrows(flow, display, 16);
        } else {
            // Lucas-Kanade sparse tracking
            if(prevPts.empty()){
                goodFeaturesToTrack(prevGray, prevPts, 1000, 0.01, 7);
            }
            if(!prevPts.empty()){
                calcOpticalFlowPyrLK(prevGray, gray, prevPts, nextPts, status, err);
                // draw
                for(size_t i=0;i<nextPts.size();++i){
                    if(!status[i]) continue;
                    Point2f p = prevPts[i];
                    Point2f q = nextPts[i];
                    line(display, p, q, Scalar(0,255,0));
                    circle(display, q, 2, Scalar(0,0,255), -1);
                }
                // prepare for next iteration
                // keep only good points
                vector<Point2f> goodPrev, goodNext;
                for(size_t i=0;i<nextPts.size();++i){
                    if(status[i]){
                        goodPrev.push_back(prevPts[i]);
                        goodNext.push_back(nextPts[i]);
                    }
                }
                prevPts.swap(goodNext);
            }
        }

        writer.write(display);

        // move to next frame
        gray.copyTo(prevGray);
        frameIdx++;
    }

    writer.release();
    cap.release();

    cout << "Done. Processed " << frameIdx << " frames. Output: " << outPath << endl;
    return 0;
}
