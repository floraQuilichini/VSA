// vsa.cpp : Defines the entry point for the console application.
//

// vsa.cpp : Defines the entry point for the console application.
//
#pragma warning (disable:4996)

#include <iostream>
#include <fstream>
#include <string>
#include <typeinfo>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh_approximation/approximate_triangle_mesh.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> Mesh;
typedef Kernel::FT FT;
typedef std::vector<std::size_t> Polygon;
typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;
typedef Mesh::Property_map<face_descriptor, std::size_t> Face_proxy_pmap;
typedef boost::property_map<Mesh, boost::vertex_point_t>::type Vertex_point_map;
namespace VSA = CGAL::Surface_mesh_approximation;
namespace PMP = CGAL::Polygon_mesh_processing;


std::string extract_ext(std::string filename)
{
	size_t i = filename.rfind('.', filename.length());
	if (i != std::string::npos) {
		return(filename.substr(i + 1, filename.length() - i));
	}

	return("");
}


bool from_file_to_polygon_soup(std::string filename, std::vector<Kernel::Point_3>& points, std::vector<Polygon>& triangular_faces, std::map<Kernel::Point_3, Kernel::Vector_3>* point_normal_map = nullptr)
{
	std::string file_ext = extract_ext(filename);
	if (file_ext.compare("ply") == 0)
	{
		// read ply file and go to the section where faces are listed
		std::ifstream file;
		file.open(filename);

		if (!file)
		{
			std::cerr << "File " << filename << " could not be opened" << std::endl;
			return false;
		}

		// find number of vertices
		std::string header_line;
		int nb_vertices = 0;
		int nb_faces = 0;
		bool flag_normal = false;

		do {
			getline(file, header_line);
			std::cout << header_line << std::endl;
			if (header_line.find("element vertex") != std::string::npos)
				nb_vertices = std::atoi((header_line.erase(0, 15)).c_str());
			else
			{
				if (header_line.find("element face") != std::string::npos)
					nb_faces = std::atoi((header_line.erase(0, 13)).c_str());
				else
				{
					if (header_line.find("property float nx") != std::string::npos)
						flag_normal = true;
				}
			}
		} while (header_line.compare("end_header") != 0);


		// get points
		Kernel::RT x, y, z, nx, ny, nz;
		int counter = 0;
		if (flag_normal)
		{
			while (counter < nb_vertices)
			{
				file >> x >> y >> z >> nx >> ny >> nz;
				points.push_back(CGAL::Point_3< Kernel >::Point_3(x, y, z));
				counter = counter + 1;
				//std::cout << "x : " << x << " y : " << y << " z : " << z << std::endl;
				if (point_normal_map)
				{
					Kernel::RT nn = sqrt(nx*nx + ny*ny + nz*nz);
					point_normal_map->insert(std::make_pair(CGAL::Point_3< Kernel >::Point_3(x, y, z), CGAL::Vector_3< Kernel >::Vector_3(nx/nn, ny/nn, nz/nn)));
				}
			}
		}
		else
		{
			while (counter < nb_vertices)
			{
				file >> x >> y >> z;
				points.push_back(CGAL::Point_3< Kernel >::Point_3(x, y, z));
				counter = counter + 1;
				//std::cout << "x : " << x << " y : " << y << " z : " << z << std::endl;
			}
		}

		// get faces
		counter = 0;
		unsigned int nb_indices_per_face, i1, i2, i3;
		while (counter < nb_faces)
		{
			file >> nb_indices_per_face >> i1 >> i2 >> i3;
			if (nb_indices_per_face != 3)
			{
				std::cout << " counter : " << counter << std::endl;
				std::cout << "nb vertices : " << nb_indices_per_face << ", indices : " << i1 << " " << i2 << " " << i3 << std::endl;
				std::cout << "non triangular polygon" << std::endl;
				return false;
			}

			std::vector<std::size_t> face{ i1, i2, i3 };
			triangular_faces.push_back(face);
			counter = counter + 1;
		}

		file.close();
		std::cout << "number of vertices " << points.size() << " theorical " << nb_vertices << " faces " << triangular_faces.size() << " theretical " << nb_faces << std::endl;
		return true;
	}

	else
	{
		// not implemented yet for other types of files
		return false;
	}

}

void write_pcd_file(std::vector<Kernel::Point_3>& points, std::string filename)
{


	int width = points.size();
	int nb_points = width;
	//Kernel::Point_3 p0 = points[0];

	std::ofstream file;
	file.precision(std::numeric_limits<float>::max_digits10);
	file.open(filename);
	// write header for point cloud where coordinates are float (on 4-bytes precision)
	file << "VERSION .7\n";
	file << "FIELDS x y z\n";
	file << "SIZE 4 4 4\n";
	file << "TYPE F F F\n";
	file << "COUNT 1 1 1\n";
	file << "WIDTH " << width << std::endl;
	file << "HEIGHT 1\n";
	file << "VIEWPOINT 0 0 0 1 0 0 0\n";
	file << "POINTS" << nb_points << std::endl;
	file << "DATA ascii\n";
	// write points
	for (int i = 0; i < nb_points; i++)
		file << points[i].x() << " " << points[i].y() << " " << points[i].z() << std::endl;
	
	
	file.close();

}



void write_ply_file(std::vector<Kernel::Point_3>& points, std::string filename, std::vector<Kernel::Vector_3>* normals_ptr = nullptr)
{


	int width = points.size();
	int nb_points = width;
	//Kernel::Point_3 p0 = points[0];

	std::ofstream file;
	file.precision(std::numeric_limits<float>::max_digits10);
	file.open(filename);
	// write header for point cloud where coordinates are float (on 4-bytes precision)
	file << "ply\n";
	file << "format ascii 1.0\n";
	file << "element vertex " << points.size() << std::endl;
	file << "property float x\n";
	file << "property float y\n";
	file << "property float z\n";
	if (normals_ptr)
	{
		file << "property float nx\n";
		file << "property float ny\n";
		file << "property float nz\n";
	}
	file << "element face 0\n";
	file << "property list uchar int vertex_indices\n";
	file << "end_header\n";
	
	// write points
	if (normals_ptr)
	{
		int i = 0;
		for (auto it = normals_ptr->begin(); it != normals_ptr->end(); it++)
		{
			file << points[i].x() << " " << points[i].y() << " " << points[i].z() << " " << (*it)[0] << " " << (*it)[1] << " " << (*it)[2] << std::endl;
			i = i + 1;
		}
	}
	else
	{
		for (int i = 0; i < nb_points; i++)
			file << points[i].x() << " " << points[i].y() << " " << points[i].z() << std::endl;
	}

	file.close();

}



int main(int argc, char **argv)
{

	if (argc < 3)
	{
		std::cout << "you must enter at least 2 arguments : the input and output filepath" << std::endl;
		return EXIT_FAILURE;
	}

	std::string input_filename = argv[1];
	std::string output_filename = argv[2];
	//double subdivision_ratio = 0.5;
	double max_nb_proxies = 200;
	if (argc == 4)
	{
		max_nb_proxies = std::stod(argv[3]);
	}

	bool write_normals = false;
	if (argc == 5)
	{
		max_nb_proxies = std::stod(argv[3]);
		if (std::strcmp(argv[4], "true") == 0 || std::strcmp(argv[4], "True") == 0)
			write_normals = true;
		else
			write_normals = false;
	}

	std::cout << "Input file:  " << argv[1] << std::endl;
	std::cout << "Output file: " << argv[2] << std::endl;

	// read input surface triangle mesh
	Mesh mesh;
	/*std::string filename = "C:/Registration/VSA/noised_source.ply";
	std::ifstream file("C:/Registration/VSA/pipe_cleaned.off");
	file >> mesh;*/

	// cast input non manifold mesh to a polygon soup
	std::vector<Kernel::Point_3> points;
	std::vector<Polygon> triangular_faces;
	std::map<Kernel::Point_3, Kernel::Vector_3> point_normal_map;
	from_file_to_polygon_soup(input_filename, points, triangular_faces, &point_normal_map);
	// repair polygon soup
	PMP::repair_polygon_soup(points, triangular_faces);
	// orient polygon soup
	PMP::orient_polygon_soup(points, triangular_faces);
	// mesh polygon soup
	PMP::polygon_soup_to_polygon_mesh(points, triangular_faces, mesh);

	std::cout << "Mesh has " << num_vertices(mesh) << " vertices and " << num_faces(mesh) << " faces" << std::endl;

	// face proxy index property map
	Face_proxy_pmap fpxmap = mesh.add_property_map<face_descriptor, std::size_t>("f:proxy_id", 0).first;

	// output planar proxies (ie normals in case we use L21 metric -- ie default case)
	std::vector<Kernel::Vector_3> proxies;

	// output indexed triangles
	std::vector<Kernel::Point_3> anchors;
	std::vector<CGAL::cpp11::array<std::size_t, 3> > triangles; // triplets of indices 
																// free function interface with named parameters
	bool is_manifold = VSA::approximate_triangle_mesh(mesh,
		CGAL::parameters::verbose_level(VSA::VERBOSE).
		seeding_method(VSA::HIERARCHICAL). // hierarchical seeding
		max_number_of_proxies(max_nb_proxies). // seeding with maximum number of proxies
											   //subdivision_ratio(subdivision_ratio). // set chord subdivision ratio threshold when meshing
		number_of_iterations(15). // number of clustering iterations after seeding
		anchors(std::back_inserter(anchors)). // anchor vertices
		triangles(std::back_inserter(triangles)). // indexed triangles
		face_proxy_map(fpxmap). // output face-proxy map
		proxies(std::back_inserter(proxies))); // output proxies
	std::cout << "#nb proxies : " << proxies.size() << std::endl;
	std::cout << "#anchor vertices: " << anchors.size() << std::endl;
	std::cout << "#triangles: " << triangles.size() << std::endl;



	// compute the barycenter of each face region
	Vertex_point_map vpmap = get(boost::vertex_point, const_cast<Mesh &>(mesh));
	std::map<size_t, std::vector<std::tuple<face_descriptor, CGAL::Point_3<CGAL::Epick>, FT>>> region_faces_barycenters;

	BOOST_FOREACH(face_descriptor f, faces(mesh))
	{
		size_t face_index = fpxmap[f];
		if (region_faces_barycenters.count(face_index) == 0) // index value has not been used yet
		{
			std::vector<std::tuple<face_descriptor, CGAL::Point_3<CGAL::Epick>, FT>> new_index_vector;   // create vector index

																										 // compute face_barycenter
			const halfedge_descriptor he = halfedge(f, mesh);
			const Kernel::Point_3 &p0 = vpmap[source(he, mesh)];
			const Kernel::Point_3 &p1 = vpmap[target(he, mesh)];
			const Kernel::Point_3 &p2 = vpmap[target(next(he, mesh), mesh)];

			CGAL::Point_3<CGAL::Epick> face_centroid = CGAL::centroid(p0, p1, p2);

			// compute face area 
			FT face_area = CGAL::sqrt(CGAL::squared_area(p0, p1, p2));

			// push into vector
			new_index_vector.push_back(std::make_tuple(f, face_centroid, face_area));

			// push into map
			region_faces_barycenters.insert(std::make_pair(face_index, new_index_vector));
		}

		else
		{
			// extract index-region vector
			std::map<size_t, std::vector<std::tuple<face_descriptor, CGAL::Point_3<CGAL::Epick>, FT>>>::iterator it = region_faces_barycenters.find(face_index);
			std::vector<std::tuple<face_descriptor, CGAL::Point_3<CGAL::Epick>, FT>> index_vector = it->second;

			// compute face_barycenter
			const halfedge_descriptor he = halfedge(f, mesh);
			const Kernel::Point_3 &p0 = vpmap[source(he, mesh)];
			const Kernel::Point_3 &p1 = vpmap[target(he, mesh)];
			const Kernel::Point_3 &p2 = vpmap[target(next(he, mesh), mesh)];

			CGAL::Point_3<CGAL::Epick> face_centroid = CGAL::centroid(p0, p1, p2);

			// compute face area 
			FT face_area = CGAL::sqrt(CGAL::squared_area(p0, p1, p2));

			// push into vector
			index_vector.push_back(std::make_tuple(f, face_centroid, face_area));

			// push into map
			region_faces_barycenters.insert(std::make_pair(face_index, index_vector));

		}

	}

/*
	// compute the weighted barycenter of each region
	std::vector<Kernel::Point_3> relevant_points;
	for (auto it = region_faces_barycenters.begin(); it != region_faces_barycenters.end(); it++)
	{
		std::vector<std::tuple<face_descriptor, CGAL::Point_3<CGAL::Epick>, FT>> index_vector = it->second;
		CGAL::Point_3<CGAL::Epick> region_barycenter = CGAL::ORIGIN;
		FT sum_areas(0.0);
		for (size_t i = 0; i < index_vector.size(); i++)
		{
			region_barycenter = region_barycenter + (std::get<1>(index_vector[i]) - CGAL::ORIGIN)*std::get<2>(index_vector[i]);
			sum_areas += std::get<2>(index_vector[i]);
		}
		Kernel::Vector_3 center = (region_barycenter - CGAL::ORIGIN) / sum_areas;
		relevant_points.push_back(CGAL::ORIGIN + center);
	}

*/
	
	// extract mesh region closest point from region barycenter
	std::vector<Kernel::Point_3> relevant_points;
	//int counter = 0;
	std::vector<Kernel::Vector_3> normals;
	for (auto it = region_faces_barycenters.begin(); it != region_faces_barycenters.end(); it++)
	{
		// compute the weighted barycenter of each region
		std::vector<std::tuple<face_descriptor, CGAL::Point_3<CGAL::Epick>, FT>> index_vector = it->second;
		CGAL::Point_3<CGAL::Epick> region_barycenter = CGAL::ORIGIN;
		FT sum_areas(0.0);
		for (size_t i = 0; i < index_vector.size(); i++)
		{
			region_barycenter = region_barycenter + (std::get<1>(index_vector[i]) - CGAL::ORIGIN)*std::get<2>(index_vector[i]);
			sum_areas += std::get<2>(index_vector[i]);
		}
		Kernel::Vector_3 center = (region_barycenter - CGAL::ORIGIN)/ sum_areas;

		// find the closest face from region center
		Kernel::FT min_dist = CGAL::squared_distance(CGAL::ORIGIN + center, std::get<1>(index_vector[0]));
		face_descriptor seed = std::get<0>(index_vector[0]);
		for (size_t i = 0; i < index_vector.size(); i++)
		{
			Kernel::FT dist = CGAL::squared_distance(CGAL::ORIGIN + center, std::get<1>(index_vector[i]));
			if (dist < min_dist)
			{
				min_dist = dist;
				seed = std::get<0>(index_vector[i]);
			}
		}

		// and then take medoid
		const halfedge_descriptor he = halfedge(seed, mesh);
		const Kernel::Point_3 &p0 = vpmap[source(he, mesh)];
		const Kernel::Point_3 &p1 = vpmap[target(he, mesh)];
		const Kernel::Point_3 &p2 = vpmap[target(next(he, mesh), mesh)];
		Kernel::Point_3 medoid;
		if (CGAL::squared_distance(CGAL::centroid(p0, p1, p2), p0) < CGAL::squared_distance(CGAL::centroid(p0, p1, p2), p1))
		{
			if (CGAL::squared_distance(CGAL::centroid(p0, p1, p2), p0) < CGAL::squared_distance(CGAL::centroid(p0, p1, p2), p2))
				medoid = p0;
			else
			{
				if (CGAL::squared_distance(CGAL::centroid(p0, p1, p2), p1) < CGAL::squared_distance(CGAL::centroid(p0, p1, p2), p2))
					medoid = p1;
				else
					medoid = p2;
			}
		}
		else
		{
			if (CGAL::squared_distance(CGAL::centroid(p0, p1, p2), p1) < CGAL::squared_distance(CGAL::centroid(p0, p1, p2), p2))
				medoid = p1;
			else
				medoid = p2;
		}
		relevant_points.push_back(medoid);
		//relevant_points.push_back(CGAL::ORIGIN + center);

		// check that <medoid-normal> are correct
		std::map<Kernel::Point_3, Kernel::Vector_3>::iterator iter = point_normal_map.find(medoid);
		normals.push_back(iter->second);
		/*std::cout << "theoretical values : " << std::endl;
		std::cout << "point : " << (iter->first).x() << " " << (iter->first).y() << " " << (iter->first).z() << " normalized normal : " << (iter->second)[0] << " " << (iter->second)[1] << " " << (iter->second)[2] << std::endl;
		std::cout << "experimental values : " << std::endl;
		std::cout << "point : " << medoid.x() << " " << medoid.y() << " " << medoid.z() << " normal : " << (proxies[counter])[0] << " " << (proxies[counter])[1] << " " << (proxies[counter])[2] << std::endl;
		counter = counter + 1;*/
	}

	

	// add anchors and middle region
	//relevant_points.insert(relevant_points.end(), anchors.begin(), anchors.end());

	// write pcd file
	//std::string pcd_filename = "C:/Registration/VSA/processed_source.pcd";
	if (write_normals)
	{
		//write_ply_file(relevant_points, output_filename, &proxies);
		write_ply_file(relevant_points, output_filename, &normals);
	}
	else
		write_pcd_file(relevant_points, output_filename);
	//write_pcd_file(anchors, output_filename);


	if (is_manifold) {
		std::cout << "oriented, 2-manifold output." << std::endl;
		// convert from soup to surface mesh
		PMP::orient_polygon_soup(anchors, triangles);
		Mesh output;
		PMP::polygon_soup_to_polygon_mesh(anchors, triangles, output);
		if (CGAL::is_closed(output) && (!PMP::is_outward_oriented(output)))
			PMP::reverse_face_orientations(output);
		/* std::ofstream out("C:/Registration/VSA/dump.off");
		out << output;
		out.close();*/
	}


	//std::cin.get();
	return EXIT_SUCCESS;
}



