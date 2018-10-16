#include <OpenGP/GL/TrackballWindow.h>
#include <OpenGP/SurfaceMesh/GL/SurfaceMeshRenderShaded.h>
#include <OpenGP/SurfaceMesh/GL/SurfaceMeshRenderFlat.h>
#include <OpenGP/GL/PointsRenderer.h>
#include <OpenGP/GL/SegmentsRenderer.h>
#include <random>

using namespace OpenGP;

struct MainWindow : public TrackballWindow {
	PointsRenderer render_points = PointsRenderer();
	SegmentsRenderer render_segments = SegmentsRenderer();

	MatMxN points;
	MatMxN points_3d_render;
	int num_points;
	double radius = 0.3;
	SegmentsRenderer::Segments segments;

	/* Used to store the original length of the polyline. Otherwise if we
	try to recompute it at every iteration of the smoothing algorithm the
	floating-point approximation error accumulates and the shapes end-up shrinking.*/
	double originalLength = 0.f;

	// ============================================================================
	// Exercise 2 : fill the 2 functions below (see PDF for instructions)
	// To test your implementation, use the S key for laplacian smoothing and the
	// C key for the osculating circle.
	// Hint : try to play with epsilon
	// ============================================================================

	// Time step for smoothing
	double epsilon = 0.1;


	/*	@brief	Computes the Euclidian distance between 2 points in 3-dimensional space
	
		@param	a	A 3D point
		@parma	b	Another 3D point

		@return	The Euclidian distance between points a and b
	*/
	double euclidianDistance(Point a, Point b) {
		return sqrt(a.x()*a.x() + b.x()*b.x());
	}

	/*	@brief	Creates a 3D point from a column of the "points" matrix

		@param	index	The index of the point to consider

		@return	A 3D point corresponding to the index-th column of the "points" matrix
	*/
	Point extendPoint(size_t index) {
		return Point(points(0, index), points(1, index), 0.f);
	}

	/*	@brief	Computes the length of the polyline represented by "points"

		@return	The length of the polyline represented by "points"
	*/
	double getPolylineLength() {

		// Iterate over the segments
		double totalLength = 0.0f;
		for (size_t i = 1; i < num_points; ++i) {
			totalLength += euclidianDistance(extendPoint(i), extendPoint(i - 1));
		}

		// Add the last segment
		totalLength += euclidianDistance(extendPoint(0), extendPoint(num_points - 1));
		return totalLength;
	}

	/*	@brief	Computes the centroid of the polyline represented by "points"

		@return	The centroid of the polyline represented by "points"
	*/
	Point getPolylineCentroid() {

		// Sum up the coordinates
		double x, y = 0.f;
		for (size_t i = 0; i < num_points; ++i) {
			x += points(0, i);
			y += points(1, i);
		}

		// Divide my the number of points and return
		return Point(x / num_points, y / num_points, 0.f);
	}

	void laplacianSmoothing() {

		// Get current centroid
		Point centroid = getPolylineCentroid();

		// Create a local copy of the polyline
		MatMxN pointsCopy = points;

		// Compute the coordinates of the new points
		points.col(0) = ((1 - epsilon) * pointsCopy.col(0)) + epsilon * ((pointsCopy.col(num_points - 1) + pointsCopy.col(1)) / 2);
		for (size_t i = 1; i < num_points - 1; ++i) {
			points.col(i) = ((1 - epsilon) * pointsCopy.col(i)) + epsilon * ((pointsCopy.col(i - 1) + pointsCopy.col(i + 1)) / 2);
		}
		points.col(num_points - 1) = ((1 - epsilon) * pointsCopy.col(num_points - 1)) + epsilon * ((pointsCopy.col(num_points - 2) + pointsCopy.col(0)) / 2);

		// Compute length ratio
		double newLength = getPolylineLength();
		double lengthRatio = originalLength / newLength;

		// Rescale the new polyline
		for (size_t i = 0; i < num_points; ++i) {

			// Compute new point
			Eigen::Vector3f vectorFromCentroid = extendPoint(i) - centroid;
			Point newPoint = centroid + (vectorFromCentroid * lengthRatio);

			// Update current point
			points(0, i) = newPoint.x();
			points(1, i) = newPoint.y();
		}

	}

	Point calculateCircumscribed(const Point &p1, const Point &p2, const Point &p3) {
		// c : circumscribed circle
		Point c = Point(0.f, 0.f, 0.f);

		float u = (p3.x()*p3.x() - p2.x()*p2.x() + p3.y()*p3.y() - p2.y()*p2.y()) / (2 * (p3.y() - p2.y()));
		float v = (p2.x()*p2.x() - p1.x()*p1.x() + p2.y()*p2.y() - p1.y()*p1.y()) / (2 * (p2.y() - p1.y()));
		float w = (((p2.x() - p1.x()) / (p2.y() - p1.y())) - ((p3.x() - p2.x()) / (p3.y() - p2.y())));

		c.x() = (u - v) / w;
		c.y() = (-(p2.x() - p1.x()) / (p2.y() - p1.y())) * c.x() + v;

		return c;
	}

	void osculatingCircle() {

		Point p1 = Point(points(0, num_points - 1), points(1, num_points - 1), 0.0f);
		Point p2 = Point(points(0, 0), points(1, 0), 0.0f);
		Point p3 = Point(points(0, 1), points(1, 1), 0.0f);

		Point c = calculateCircumscribed(p1, p2, p3);

		Point newPoint = p2 + epsilon * (c - p2) / (c - p2).squaredNorm();

		points(0, 0) = newPoint.x();
		points(1, 0) = newPoint.y();

		for (int i = 1; i < num_points - 1; ++i) {
			p1 = Point(points(0, i - 1), points(1, i - 1), 0.0f);
			p2 = Point(points(0, i), points(1, i), 0.0f);
			p3 = Point(points(0, i + 1), points(1, i + 1), 0.0f);

			c = calculateCircumscribed(p1, p2, p3);

			newPoint = p2 + epsilon * (c - p2) / (c - p2).squaredNorm();

			points(0, i) = newPoint.x();
			points(1, i) = newPoint.y();
		}

		p1 = Point(points(0, num_points - 2), points(1, num_points - 2), 0.0f);
		p2 = Point(points(0, num_points - 1), points(1, num_points - 1), 0.0f);
		p3 = Point(points(0, 0), points(1, 0), 0.0f);

		c = calculateCircumscribed(p1, p2, p3);

		newPoint = p2 + epsilon * (c - p2) / (c - p2).squaredNorm();

		points(0, num_points - 1) = newPoint.x();
		points(1, num_points - 1) = newPoint.y();

		// std::cout << c.x() << " " << c.y() << " " << c.z() << std::endl;


		// RESCALE HERE

		// Curve Smoothing - osculating circle (again, this function should do one iteration of smoothing)
	}

	// ============================================================================
	// END OF Exercise 2 (do not thouch the rest of the code)
	// ============================================================================

	void generateRandomizedClosedPolyline() {
		std::default_random_engine generator;
		std::uniform_real_distribution<double> distribution(0., 5 * 3e-2);

		Vec2 center(3e-2, 2e-3);

		points = MatMxN::Zero(2, num_points);
		for (int i = 0; i < num_points; ++i)
		{
			double frac = static_cast<double>(i) / static_cast<double>(num_points);
			points(0, i) = center(0) + radius * cos(2. * M_PI * frac) + distribution(generator);
			points(1, i) = center(1) + radius * sin(2. * M_PI * frac) + distribution(generator);
		}

		/* We remember the original length of the polyline here. Otherwise if we
		try to recompute it at every iteration of the smoothing algorithm the 
		floating-point approximation error accumulates and the shapes end-up shrinking.*/
		originalLength = getPolylineLength();

	}

	void render() {

		// Prepare the render points
		points_3d_render = MatMxN::Zero(3, points.cols());
		points_3d_render.block(0, 0, 2, points.cols()) = points;

		// Rebuild the segments
		segments.clear();
		for (int i = 0; i < points_3d_render.cols(); ++i) {
			segments.push_back({ points_3d_render.col(i), points_3d_render.col((i + 1) % points_3d_render.cols()) });
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
			// TODO: remove!
			for (int i = 0; i < 10; ++i) {
				laplacianSmoothing();
			}
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
