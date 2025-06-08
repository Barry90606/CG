#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
using namespace std;

const int width = 1000;
const int height = 1000;

struct Color {
    int r, g, b;
    Color(int r = 255, int g = 255, int b = 255): r(r), g(g), b(b) {}
};

struct Point {
    double x, y;
    Point(double x = 0, double y = 0): x(x), y(y) {}
    Point operator+(const Point& p) const { return Point(x + p.x, y + p.y); }
    Point operator-(const Point& p) const { return Point(x - p.x, y - p.y); }
    Point operator*(double s) const { return Point(x * s, y * s); }
    Point operator/(double s) const { return Point(x / s, y / s); }
};

// Calculate dot product
inline double dot(const Point& a, const Point& b) {
    return a.x * b.x + a.y * b.y;
}

// Calculate distance between two points
inline double distance(const Point& a, const Point& b) {
    double dx = a.x - b.x;
    double dy = a.y - b.y;
    return std::sqrt(dx * dx + dy * dy);
}

// Compute area of triangle (Heron's formula)
double triangleArea(const Point& a, const Point& b, const Point& c) {
    double ab = distance(a, b);
    double bc = distance(b, c);
    double ca = distance(c, a);
    double s = 0.5 * (ab + bc + ca);
    return sqrt(s * (s - ab) * (s - bc) * (s - ca));
}

double computenewlambda(const Point& ci0, const Point& ci1, const Point& ci11, const Point& ci12){
    double area1 = sqrt(triangleArea(ci0, ci1, ci11));
    double area2 = sqrt(triangleArea(ci0, ci1, ci11));
    double area3 = sqrt(triangleArea(ci1, ci11, ci12));
    double lambda = area1 / (area2 + area3);
    return lambda;

    // return 0.5;
}

// Solve cubic function root using binary search
double solveT(const Point& ci0, const Point& ci2, const Point& pi) {
    // (ci2 - ci0) * (ci2 - ci0) * t^3 + 3 * (ci2 - ci0) * (ci0 - pi) * t^2 + (3 * ci0 - 2 * pi - ci2) * (ci0 -pi) * t - (ci0 - pi) * (ci0 - pi) = 0
    const double delta = 1e-4;
    Point c20 = ci2 - ci0;
    Point c0_pi = ci0 - pi;

    double A = dot(c20, c20);
    double B = 3.0 * dot(c0_pi, c20);

    Point p1 = ci0 * 3.0 - pi * 2.0 - ci2;
    Point p2 = ci0 - pi;

    double C = dot(p1, p2);
    double D = -dot(p2, p2);

    double t0 = 0.0, t1 = 1.0, t;
    while (t1 - t0 > delta) {
        t = 0.5 * (t0 + t1);
        // cout << t <<endl;
        double f = A * t * t * t + B * t * t + C * t + D;
        double f0 = A * t0 * t0 * t0 + B * t0 * t0 + C * t0 + D;
        //如果符號相反，表示根在[t0​ ,t]，否則在[t,t1​ ]
        //意思是 f(a) 和 f(b) 一正一負，那麼根據介值定理，這段區間內必定至少存在一個根 t∈(a,b)，使得 f(t)=0。
        if(f == 0){
            cout << "time1" <<endl;
            return t;
        }
        if(f0 == 0){
            cout << "time2" <<endl;
            return t0;
        }
        if (f * f0 < 0) t1 = t;
        else t0 = t;    
    }
    cout << "timeout" <<endl;
    return 0.5 * (t0 + t1);
}

// Compute middle control point Ci1
Point computeCi1(const Point& pi, double t, const Point& ci0, const Point& ci2) {
    double denom = 2 * t * (1 - t);
    if (fabs(denom) < 1e-6) denom = 1e-6;
    Point ci1 = (pi - ci0 * (1 - t) * (1 - t) - ci2 * t * t) / denom;
    // Point ci1 = (pi - ci0 * (1 - t) * (1 - t) - ci2 * t * t) / (2 * t * (1 - t));
    return ci1;
}

Point computeCi2(double lambda, const Point& ci1, const Point& ci11){
    Point ci2 = ci1 * (1 - lambda) + ci11 * lambda;
    return ci2;
}

Color image[height][width];

void draw_point(int x, int y, Color c = Color(0, 0, 0), int size = 1) {
    int half = size / 2;
    for (int dx = -half; dx <= half; ++dx)
        for (int dy = -half; dy <= half; ++dy) {
            int nx = x + dx, ny = y + dy;
            if (nx >= 0 && nx < width && ny >= 0 && ny < height)
                image[height - 1 - ny][nx] = c;
        }
}

void draw_bezier(Point p0, Point p1, Point p2, Color c = Color(0, 0, 0)) {
    for (double t = 0; t <= 1.0; t += 0.001) {
        double x = (1 - t)*(1 - t)*p0.x + 2*(1 - t)*t*p1.x + t*t*p2.x;
        double y = (1 - t)*(1 - t)*p0.y + 2*(1 - t)*t*p1.y + t*t*p2.y;
        draw_point(round(x), round(y), c);
    }
}

void clear_image() {
    for (int j = 0; j < height; ++j)
        for (int i = 0; i < width; ++i)
            image[j][i] = Color(255, 255, 255);
}

void save_image(const string& filename) {
    ofstream file(filename);
    file << "P3\n" << width << " " << height << "\n255\n";
    for (int j = 0; j < height; ++j)
        for (int i = 0; i < width; ++i)
            file << image[j][i].r << ' ' << image[j][i].g << ' ' << image[j][i].b << '\n';
    file.close();
}

int main() {
    clear_image();

    vector<Point> A = {
        Point(60, 60),
        Point(140, 60),
        Point(100, 140)
    };


    int N = A.size();
    cout << N << endl;

    vector<double> lambda;
    for( int l = 0; l < N; l++){
        lambda.push_back(0.5);
    }
    // 畫資料點
    for (const auto& p : A)
        draw_point(round(p.x), round(p.y), Color(0, 0, 255), 3);

    vector<Point> inputs_ci1;
    for( int l = 0; l < N; l++){
        inputs_ci1.push_back(A[l]);
        // cout<<inputs_ci1[l].x<<" "<<inputs_ci1[l].y<<endl;
    }
    vector<Point> inputs_ci2;
    for(int i = 0; i < N ;i++){
        Point ci1 = inputs_ci1[i];
        Point ci11 = inputs_ci1[(i+1)%N]; 
        double temp_lambda = lambda[i];
        Point ci2 = computeCi2(temp_lambda, ci1, ci11);
        inputs_ci2.push_back(ci2);
    }
    // for(int i = 0; i < N; i++){
    //     cout<<i<<endl;
    //     cout<<A[i].x<<" "<<A[i].y<<endl;
    //     cout<<inputs_ci1[i].x<<" "<<inputs_ci1[i].y<<endl;
    //     cout<<inputs_ci2[i].x<<" "<<inputs_ci2[i].y<<endl;
    // }
    int time = 5;
    for( int t = 0; t < time ; t++){
        // ci1
        for( int j = 0; j < N ; j++){
            Point ci0 = inputs_ci2[(j+N-1)%N];
            Point ci2 = inputs_ci2[j];
            double t = solveT(ci0, ci2, A[j]);
            Point new_ci1 = computeCi1(A[j], t, ci0, ci2);
            inputs_ci1[j] = new_ci1;
        }
        // lambda
        for( int k = 0; k < N; k++){
            Point ci0 = inputs_ci2[(k+N-1)%N];
            Point ci1 = inputs_ci1[k];
            Point ci11 = inputs_ci1[(k+1)%N];
            Point ci12 = inputs_ci2[(k+1)%N];
            lambda[k] = computenewlambda(ci0, ci1, ci11, ci12);
        }
        // ci2
        for( int h = 0; h < N; h++){
            Point ci1 = inputs_ci1[h];
            Point ci11 = inputs_ci1[(h+1)%N];
            double temp_lambda = lambda[h];
            Point new_ci2 = computeCi2(temp_lambda, ci1, ci11);
            inputs_ci2[h] = new_ci2;
        }
    }

    for (int i = 0; i < N; i++) {
        cout<<i<<endl;
        Point ci0 = inputs_ci2[(i-1+N)%N];
        Point ci1 = inputs_ci1[i];
        Point ci2 = inputs_ci2[i];
        cout<<ci0.x<<"  "<<ci0.y<<endl;
        cout<<ci1.x<<"  "<<ci1.y<<endl;
        cout<<ci2.x<<"  "<<ci2.y<<endl;
        draw_point(round(ci0.x), round(ci0.y), Color(0, 255, 255), 3);
        draw_point(round(ci1.x), round(ci1.y), Color(0, 255, 255), 3);
        draw_point(round(ci2.x), round(ci2.y), Color(0, 255, 255), 3);
        draw_bezier(ci0, ci1, ci2);
    }

    save_image("test2_8.ppm");
    cout << "Saved image to closed_bezier_correct.ppm\n";
    return 0;
}
