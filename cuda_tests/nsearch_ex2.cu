#include "../lib/cuNSearch/src/cuNSearchDeviceData.h"

using namespace cuNSearch;

//using Real3 = std::array<Real, 3>;
std::size_t const N = 50;
//std::vector<Real3> positions;
std::vector< std::array<Real, 3> > positions;
Real const r_omega = static_cast<Real>(0.01/2.)/ static_cast<Real>(N - 1);
Real const radius = static_cast<Real>(2.) * (static_cast<Real>(2.0) * r_omega );

double3 *x_d; 

int main() {
	
	//cuNSearch::NeighborhoodSearch nsearch(r);

	//Generate test data
	Real min_x = std::numeric_limits<Real>::max();
	Real max_x = std::numeric_limits<Real>::min();
	positions.reserve(N * N * N);
	for (unsigned int i = 0; i < N; ++i)
	{
		for (unsigned int j = 0; j < N; ++j) {
			for (unsigned int k = 0; k < N; ++k) {
				std::array<Real, 3> x = { {
						// r_omega * static_cast<Real>(2.0 * static_cast<double>(i) / static_cast<double>(N - 1) - 1.0),
						// r_omega * static_cast<Real>(2.0 * static_cast<double>(j) / static_cast<double>(N - 1) - 1.0),
						// r_omega * static_cast<Real>(2.0 * static_cast<double>(k) / static_cast<double>(N - 1) - 1.0) 
						
						r_omega * static_cast<Real>(2.0 * static_cast<double>(i) - 1.0),
						r_omega * static_cast<Real>(2.0 * static_cast<double>(j) - 1.0),
						r_omega * static_cast<Real>(2.0 * static_cast<double>(k) - 1.0) 						
						
						
						} };
				
				positions.push_back(x);
			}
		}
	}

	std::random_shuffle(positions.begin(), positions.end());
	printf("Number of particles: %d \n", static_cast<int>(positions.size()));

	//Create neighborhood search instance
	NeighborhoodSearch nsearch(radius);

	//Add point set from the test data
	auto pointSetIndex = nsearch.add_point_set(positions.front().data(), positions.size(), true, true);

	for (size_t i = 0; i < 5; i++)
	{
		if (i != 0)
		{
			nsearch.z_sort();
			nsearch.point_set(pointSetIndex).sort_field((Real3*)nsearch.point_set(pointSetIndex).GetPoints());
		}

		Timing::reset();
		nsearch.find_neighbors();
		Timing::printAverageTimes();
	}

	//Neighborhood search result test
	auto &pointSet = nsearch.point_set(0);
	auto points = pointSet.GetPoints();
	
	//Obtaining and update device particle positions
	// auto &pointSet = nsearch.point_set(0);
	// auto points = pointSet.GetPoints();
	
	cout << "Allocating" << endl;
	cudaMalloc((void **)&x_d, N * N * N * sizeof (double3));
	cudaMemcpy(x_d, x, size, cudaMemcpyHostToDevice);
	
	
	// unsigned int point_set_id = nsearch.add_point_set(positions.front().data(), positions.size());
	//nsearch.find_neighbors();	
	//
	//PointSet(x, n, is_dynamic, user_data default null)
	//PointSet ps(positions.front().data(), positions.size(), true);
	
	
	//cuNSearchDeviceData data_d;
	
}
	//Accesing points (from demo)

	// //Neighborhood search result test
	// auto &pointSet = nsearch.point_set(0);
	// auto points = pointSet.GetPoints();

	// std::cout << "Validate results" << std::endl;
	// for (unsigned int i = 0; i < pointSet.n_points(); i++)
	// {
		// Real3 point = ((Real3*)points)[i];
		// auto count = pointSet.n_neighbors(0, i);
		// for (unsigned int j = 0; j < count; j++)
		// {
			// auto neighbor = pointSet.neighbor(0, i, j);
			// auto diff = point - ((Real3*)points)[neighbor];
			// float squaredLength = diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2];
			// float distance = sqrt(squaredLength);

			// if (distance > radius)
			// {
				// throw std::runtime_error("Not a neighbor");
			// }
		// }
	// }


	
		// unsigned int NeighborhoodSearch::add_point_set(Real const* x, std::size_t n, bool is_dynamic,
		// bool search_neighbors, bool find_neighbors, void *user_data)
	// {
		// auto index = pointSets.size();
		// pointSets.push_back(PointSet(x, n, is_dynamic, user_data));
		// m_activation_table.add_point_set(search_neighbors, find_neighbors);

		// for (auto &pointSet : pointSets)
		// {
			// pointSet.neighbors.resize(pointSets.size());
		// }

		// return static_cast<unsigned int>(index);
	// }

	// void
		// NeighborhoodSearch::find_neighbors(bool points_changed_)
	// {
		// if (points_changed_ || !isInitialized)
		// {
			// for (auto &pointSet : pointSets)
			// {
				// if (!isInitialized || pointSet.is_dynamic())
				// {
					// updatePointSet(pointSet);
				// }
			// }
		// }
		// isInitialized = true;

		// for (unsigned int i = 0; i < pointSets.size(); i++)
		// {
			// for (unsigned int j = 0; j < pointSets.size(); j++)
			// {
				// if (m_activation_table.is_active(i, j))
				// {
					// auto &queryPointSet = pointSets[i];
					// auto &pointSet = pointSets[j];
					// deviceData->computeNeighborhood(queryPointSet, pointSet, j);
				// }
			// }
		// }
	// }