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

	// Time step for Laplacian smoothing
	double epsilonLaplacian = 0.1;

	// Time step for osculating circle smoothing
	double epsilonOsculating = 0.00005;

	/*	@brief	Computes the Euclidian distance between 2 points in 3-dimensional space

		@param	a	A 3D point
		@parma	b	Another 3D point

		@return	The Euclidian distance between points a and b
	*/
	double euclidianDistance(Point a, Point b) {
		return sqrt(a.x()*a.x() + b.x()*b.x());
	}

	/*	@brief	Creates a 3D point from a column of the polyline matrix

		@param	polyline	The polyline to consider
		@param	index		The index of the point to consider

		@return	A 3D point corresponding to the index-th column of the polyline matrix
	*/
	Point extendPoint(const MatMxN &polyline, size_t index) {
		return Point(polyline(0, index), polyline(1, index), 0.f);
	}

	/*	@brief	Creates a 3D point from a column of the "points" matrix

		@param	index		The index of the point to consider

		@return	A 3D point corresponding to the index-th column of the "points" matrix
	*/
	Point extendPoint(size_t index) {
		return extendPoint(points, index);
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
		float x = 0.f, y = 0.f;
		for (size_t i = 0; i < num_points; ++i) {
			x += points(0, i);
			y += points(1, i);
		}

		// Divide my the number of points and return
		return Point(x / num_points, y / num_points, 0.f);
	}

	/*	@brief	Rescales the polyline "points" by a certain factor

		@param	scalingFactor	The scaling factor
		@return	void
	*/
	void rescalePolyline(double scalingFactor) {

		// Get current centroid
		Point centroid = getPolylineCentroid();

		// Rescale the new polyline
		for (size_t i = 0; i < num_points; ++i) {

			// Compute new point
			Eigen::Vector3f vectorFromCentroid = extendPoint(i) - centroid;
			Point newPoint = centroid + (vectorFromCentroid * scalingFactor);

			// Update current point
			points(0, i) = newPoint.x();
			points(1, i) = newPoint.y();
		}
	}

	/*	@brief	Updates the coordinates of a point using Laplacian smoothing

		@param	polyline	The polyline of reference
		@param	indexLow	Index of previous vertex
		@param	indexCur	Index of current vertex (the one being updated)
		@param	indexHigh	Index of next vertex
	*/
	void updatePointLaplacian(const MatMxN &polyline, size_t indexLow, size_t indexCur, size_t indexHigh) {
		points.col(indexCur) = ((1 - epsilonLaplacian) * polyline.col(indexCur)) +
			epsilonLaplacian * ((polyline.col(indexLow) + polyline.col(indexHigh)) / 2);
	}

	void laplacianSmoothing() {

		// Create a local copy of the polyline
		MatMxN pointsCopy = points;

		// Compute the coordinates of the new points
		updatePointLaplacian(pointsCopy, num_points - 1, 0, 1);
		for (size_t i = 1; i < num_points - 1; ++i) {
			updatePointLaplacian(pointsCopy, i - 1, i, i + 1);
		}
		updatePointLaplacian(pointsCopy, num_points - 2, num_points - 1, 0);

		// Rescale the polyline
		rescalePolyline(originalLength / getPolylineLength());
	}

	/*	@brief	Computes the center coordinates of the circumscrubed circle to triangle ABC

		@param	a	Point A
		@param	b	Point B
		@param	c	Point C

		@return	The center coordinates of the circumscrubed circle to triangle ABC
	*/
	Point computeCircumscribed(const Point &a, const Point &b, const Point &c) {

		double d = 2 * (
			a.x() * (b.y() - c.y()) +
			b.x() * (c.y() - a.y()) +
			c.x() * (a.y() - b.y())
			);
		double u = (1.f / d) * (
			(a.x()*a.x() + a.y()*a.y()) * (b.y() - c.y()) +
			(b.x()*b.x() + b.y()*b.y()) * (c.y() - a.y()) +
			(c.x()*c.x() + c.y()*c.y()) * (a.y() - b.y())
			);
		double v = (1.f / d) * (
			(a.x()*a.x() + a.y()*a.y()) * (c.x() - b.x()) +
			(b.x()*b.x() + b.y()*b.y()) * (a.x() - c.x()) +
			(c.x()*c.x() + c.y()*c.y()) * (b.x() - a.x())
			);

		return Point(u, v, 0.f);
	}

	/*	@brief	Updates the coordinates of a point using oscularing circle smoothing

		@param	polyline	The polyline of reference
		@param	indexLow	Index of previous vertex
		@param	indexCur	Index of current vertex (the one being updated)
		@param	indexHigh	Index of next vertex

		@return	void
	*/
	void updatePointOsculating(const MatMxN &polyline, size_t indexLow, size_t indexCur, size_t indexHigh) {

		Point v = extendPoint(polyline, indexCur);
		Point center = computeCircumscribed(
			extendPoint(polyline, indexLow),
			v,
			extendPoint(polyline, indexHigh)
		);

		Point newPoint = v + (epsilonOsculating * ((center - v) / (center - v).squaredNorm()));
		points(0, indexCur) = newPoint.x();
		points(1, indexCur) = newPoint.y();
	}

	void osculatingCircle() {

		// Create a local copy of the polyline
		MatMxN pointsCopy = points;

		double length = getPolylineLength();

		// Compute the coordinates of the new points
		updatePointOsculating(pointsCopy, num_points - 1, 0, 1);
		for (size_t i = 1; i < num_points - 1; ++i) {
			updatePointOsculating(pointsCopy, i - 1, i, i + 1);
		}
		updatePointOsculating(pointsCopy, num_points - 2, num_points - 1, 0);

		// Rescale the polyline
		double scalingFactor = originalLength / getPolylineLength();
		if (scalingFactor > 1) { // Only scale up
			rescalePolyline(scalingFactor);
		}
	}

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
			laplacianSmoothing();
		}
		else if (key == GLFW_KEY_C && action == GLFW_RELEASE)
		{
			std::cout << "Length: " << getPolylineLength() << std::endl;
			for (size_t i = 0 ; i < 5 ; ++i)
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
