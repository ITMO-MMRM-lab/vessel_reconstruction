Данный модуль принимает:
- файл с объемной (? или поверхностной ?) сеткой, 
- список размеченных точек (0 - фиксированная, 1 - свободная, 2 - управляющая), 
- список индекс управляющей точки - новые координаты XYZ.
[возможно целееобразней размечать точки уже в cgal? но это не точно =з]

На выходе получеам деформированную поверхность.

Для данного модуля требуются:
- gcc,
- Qt5,
- CGAL,
- OpenMesh,
- Eigen,
- Boost.

Для запуска:
..vessel_reconstruction/stent_deformation/$ mkdir build
..vessel_reconstruction/stent_deformation/$ cd build
..vessel_reconstruction/stent_deformation/build$ cmake .
..vessel_reconstruction/stent_deformation/build$ make
//рекомендую использовать cmake GUI для уточнения зависимостей