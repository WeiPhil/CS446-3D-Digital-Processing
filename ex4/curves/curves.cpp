#include <OpenGP/GL/TrackballWindow.h>
#include <OpenGP/SurfaceMesh/GL/SurfaceMeshRenderShaded.h>
#include <OpenGP/SurfaceMesh/GL/SurfaceMeshRenderFlat.h>
#include <OpenGP/GL/PointsRenderer.h>
#include <OpenGP/GL/SegmentsRenderer.h>
#include <random>

using namespace OpenGP;

struct MainWindow : public TrackballWindow {
  PointsRenderer render_points = PointsRenderer ();
  SegmentsRenderer render_segments = SegmentsRenderer ();

  MatMxN points;
  MatMxN points_3d_render;
  int num_points;
  double radius = 0.3;
  SegmentsRenderer::Segments segments;

  // !!!!!! ADDED AND CHANGED IN GIVEN CODE !!!!
  double lengthBeforeSmoothing = 0.0f;

// ============================================================================
// Exercise 2 : fill the 2 functions below (see PDF for instructions)
// To test your implementation, use the S key for laplacian smoothing and the
// C key for the osculating circle.
// Hint : try to play with epsilon
// ============================================================================
  // time step for smoothing
  double epsilon = 0.1;

  double euclidianDistance(Point a, Point b){
    return sqrt(a.x()*a.x() + b.x()*b.x());
  }

  double calculateLength(const MatMxN &points){

    Point initialPoint = Point(points(0,0), points(1,0),0.f);

    double totalLength = 0.0f;
    for (int i = 1; i < num_points; ++i){
      Point a = Point(points(0,i-1), points(1,i-1) ,0.f);
      Point b = Point(points(0,i), points(1,i),0.f );
      totalLength += euclidianDistance(a,b);
    }
    totalLength += euclidianDistance(initialPoint ,Point( points(0,num_points-1),points(1,num_points-1),0.f));

    return totalLength;
  }

  void laplacianSmoothing() {

    MatMxN pointsCopy = points;

    

    Point centroid = Point(0.0f,0.0f,0.0f);
    for (int i = 0; i < num_points; ++i){
      centroid.x() += points(0,i);
      centroid.y() += points(1,i);
    }
    centroid.x() /= num_points;
    centroid.y() /= num_points;

    std::cout << "Centroid at : " << centroid.x() << " " << centroid.y() << std::endl;

    // Calculating 
    points.col(0) = (1-epsilon) * pointsCopy.col(0) + epsilon * ( pointsCopy.col(num_points-1) + pointsCopy.col(1) )/2 ;
    for (int i = 1; i < num_points-1; ++i){
      points.col(i) = (1-epsilon) * pointsCopy.col(i) + epsilon * ( pointsCopy.col(i-1) + pointsCopy.col(i+1) )/2 ;
    }
    points.col(num_points-1) = (1-epsilon) * pointsCopy.col(num_points-1) + epsilon * ( pointsCopy.col(num_points-2) + pointsCopy.col(0) )/2 ;

    double lengthAfterSmoothing = calculateLength(points);

    double ratio = lengthBeforeSmoothing/lengthAfterSmoothing ;

    

    for (int i = 0; i < num_points; ++i){
      Eigen::Vector3f vectorFromCentroid = Point(points(0,i),points(1,i),0.0f) - centroid ;

      Point newPoint = centroid + vectorFromCentroid*ratio ;
      points(0,i) = newPoint.x();
      points(1,i) = newPoint.y();
    }
    
    std::cout << "Length of Curve : " << lengthAfterSmoothing << std::endl;

    // Curve Smoothing - centroid (this function should do one iteration of smoothing)
  }

  Point calculateCircumscribed(const Point &p1,const Point &p2,const Point &p3){
    // c : circumscribed circle
    Point c = Point(0.f,0.f,0.f);

    float u =  ( p3.x()*p3.x() - p2.x()*p2.x() + p3.y()*p3.y() - p2.y()*p2.y()) / (2*(p3.y()-p2.y()) ) ;
    float v =  ( p2.x()*p2.x() - p1.x()*p1.x() + p2.y()*p2.y() - p1.y()*p1.y()) / (2*(p2.y()-p1.y()) ) ;
    float w = ( ((p2.x()-p1.x())/(p2.y()-p1.y())) - ((p3.x() - p2.x())/(p3.y()-p2.y())) );
    
    c.x() = ( u - v ) / w ;
    c.y() = ( -(p2.x()-p1.x())/(p2.y()-p1.y()) ) * c.x() + v ;

    return c;
  }

  void osculatingCircle() {

    Point p1 = Point(points(0,num_points-1),points(1,num_points-1),0.0f);
    Point p2 = Point(points(0,0),points(1,0),0.0f);
    Point p3 = Point(points(0,1),points(1,1),0.0f);

    Point c = calculateCircumscribed(p1,p2,p3);

    Point newPoint = p2 + epsilon * (c-p2)/(c-p2).squaredNorm();

    points(0,0) = newPoint.x();
    points(1,0) = newPoint.y();

    for (int i = 1; i < num_points-1 ; ++i){
      p1 = Point(points(0,i-1),points(1,i-1),0.0f);
      p2 = Point(points(0,i),points(1,i),0.0f);
      p3 = Point(points(0,i+1),points(1,i+1),0.0f);

      c = calculateCircumscribed(p1,p2,p3);

      newPoint = p2 + epsilon * (c-p2)/(c-p2).squaredNorm();

      points(0,i) = newPoint.x();
      points(1,i) = newPoint.y();
    }

    p1 = Point(points(0,num_points-2),points(1,num_points-2),0.0f);
    p2 = Point(points(0,num_points-1),points(1,num_points-1),0.0f);
    p3 = Point(points(0,0),points(1,0),0.0f);

    c = calculateCircumscribed(p1,p2,p3);

    newPoint = p2 + epsilon * (c-p2)/(c-p2).squaredNorm();

    points(0,num_points-1) = newPoint.x();
    points(1,num_points-1) = newPoint.y();

    // std::cout << c.x() << " " << c.y() << " " << c.z() << std::endl;


    // RESCALE HERE

    // Curve Smoothing - osculating circle (again, this function should do one iteration of smoothing)
  }

// ============================================================================
// END OF Exercise 2 (do not thouch the rest of the code)
// ============================================================================

  void generateRandomizedClosedPolyline() {
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0., 5*3e-2);

    Vec2 center(3e-2, 2e-3);


    points = MatMxN::Zero(2, num_points);
    for (int i = 0; i < num_points; ++i)
    {
      double frac = static_cast<double>(i) / static_cast<double>(num_points);
      points(0, i) = center(0) + radius * cos (2. * M_PI * frac) + distribution(generator);
      points(1, i) = center(1) + radius * sin (2. * M_PI * frac) + distribution(generator);
    }

    lengthBeforeSmoothing = calculateLength(points);

  }

  void render () {

    // Prepare the render points
    points_3d_render = MatMxN::Zero(3, points.cols());
    points_3d_render.block(0, 0, 2, points.cols()) = points;

    // Rebuild the segments
    segments.clear();
    for (int i = 0; i < points_3d_render.cols(); ++i) {
      segments.push_back({ points_3d_render.col(i), points_3d_render.col((i+1) % points_3d_render.cols()) });
    }
    render_points.init_data(points_3d_render);
    render_segments.init_data(segments);
  }

  MainWindow(int argc, char** argv) : TrackballWindow("2D Viewer", 640, 480) {
    num_points = 30;
    generateRandomizedClosedPolyline();

    this->scene.add(render_points);
    this->scene.add(render_segments);

    render();
  }

  bool key_callback(int key, int scancode, int action, int mods) override {
    TrackballWindow::key_callback(key, scancode, action, mods);
    if (key == GLFW_KEY_S && action == GLFW_RELEASE)
    {
      laplacianSmoothing();     
    }
    else if (key == GLFW_KEY_C && action == GLFW_RELEASE)
    {
      osculatingCircle();
    }
    else if (key == GLFW_KEY_1 && action == GLFW_RELEASE)
    {
      num_points = 30;
      radius = 0.3;
      generateRandomizedClosedPolyline();
    }
    else if (key == GLFW_KEY_2 && action == GLFW_RELEASE)
    {
      num_points = 50;
      radius = 0.1;
      generateRandomizedClosedPolyline();
    }
    else if (key == GLFW_KEY_3 && action == GLFW_RELEASE)
    {
      num_points = 100;
      radius = 0.1;
      generateRandomizedClosedPolyline();
    }
    else if (key == GLFW_KEY_4 && action == GLFW_RELEASE)
    {
      num_points = 150;
      radius = 0.1;
      generateRandomizedClosedPolyline();
    }

    render();
    return true;
  }
};


int main(int argc, char** argv)
{
  MainWindow window(argc, argv);
  return window.run();
}
