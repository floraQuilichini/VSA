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
typedef std::vector<std::size_t> Polygon;
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


bool from_file_to_polygon_soup(std::string filename, std::vector<Kernel::Point_3>& points, std::vector<Polygon>& triangular_faces)
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

		do {
			getline(file, header_line);
			//std::cout << header_line << std::endl;
			if (header_line.find("element vertex") != std::string::npos)
				nb_vertices = std::atoi((header_line.erase(0, 15)).c_str());
			else 
			{
				if (header_line.find("element face") != std::string::npos)
					nb_faces = std::atoi((header_line.erase(0, 13)).c_str());
			}
		} while (header_line.compare("end_header") != 0);
		

		// get points
		Kernel::RT x, y, z;
		int counter = 0;
		while (counter < nb_vertices)
		{
			file >> x >> y >> z;
			points.push_back(CGAL::Point_3< Kernel >::Point_3(x, y, z));
			counter = counter + 1;
			//std::cout << "x : " << x << " y : " << y << " z : " << z << std::endl;
		}

		// get faces
		counter = 0;
		unsigned int nb_indices_per_face, i1, i2, i3;
		while (counter < nb_faces)
		{
			file >> nb_indices_per_face >> i1 >> i2 >> i3;
			if (nb_indices_per_face != 3)
			{
				std::cout << "non triangular polygon" << std::endl;
				return false;
			}

			std::vector<std::size_t> face{i1, i2, i3};
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
	{
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
	double subdivision_ratio = 0.5;
	if (argc == 4)
	{
		subdivision_ratio = std::stod(argv[3]);
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
	from_file_to_polygon_soup(input_filename, points,  triangular_faces);
		// repair polygon soup
	PMP::repair_polygon_soup(points, triangular_faces);
		// orient polygon soup
	PMP::orient_polygon_soup(points, triangular_faces);
		// mesh polygon soup
	PMP::polygon_soup_to_polygon_mesh(points, triangular_faces, mesh);

	std::cout << "Mesh has " << num_vertices(mesh) << " vertices and " << num_faces(mesh) << " faces" << std::endl;


	// output indexed triangles
	std::vector<Kernel::Point_3> anchors;
	std::vector<CGAL::cpp11::array<std::size_t, 3> > triangles; // triplets of indices
																// free function interface with named parameters
	bool is_manifold = VSA::approximate_triangle_mesh(mesh,
		CGAL::parameters::verbose_level(VSA::VERBOSE).
		seeding_method(VSA::HIERARCHICAL). // hierarchical seeding
		max_number_of_proxies(200). // seeding with maximum number of proxies
		subdivision_ratio(subdivision_ratio). // set chord subdivision ratio threshold when meshing
		number_of_iterations(30). // number of clustering iterations after seeding
		anchors(std::back_inserter(anchors)). // anchor vertices
		triangles(std::back_inserter(triangles))); // indexed triangles
	std::cout << "#anchor vertices: " << anchors.size() << std::endl;
	std::cout << "#triangles: " << triangles.size() << std::endl;


	// write pcd file
	//std::string pcd_filename = "C:/Registration/VSA/processed_source.pcd";
	write_pcd_file(anchors, output_filename);




	/*if (is_manifold) {
		std::cout << "oriented, 2-manifold output." << std::endl;
		// convert from soup to surface mesh
		PMP::orient_polygon_soup(anchors, triangles);
		Mesh output;
		PMP::polygon_soup_to_polygon_mesh(anchors, triangles, output);
		if (CGAL::is_closed(output) && (!PMP::is_outward_oriented(output)))
			PMP::reverse_face_orientations(output);
		std::ofstream out("C:/Registration/VSA/dump.off");
		out << output;
		out.close();
	}*/


//	std::cin.get();
	return EXIT_SUCCESS;
}



