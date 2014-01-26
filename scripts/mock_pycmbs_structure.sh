#!/bin/bash

function show_tree_output {
    tree $1 --noreport --dirsfirst | sed "s/^\/tmp.*/pyCMBS-v1.00/g" | awk '{ print "    "$0 }' 
}

function side_by_side_layout {

    echo ""
    echo "------"
    echo ""
    root_dir=/tmp/mock/pyCMBS-vX.YY

    mkdir $root_dir -p
    touch $root_dir/MANIFEST $root_dir/setup.py $root_dir/Makefile $root_dir/INSTALL \
            $root_dir/README $root_dir/LICENSE

    mkdir $root_dir/docs $root_dir/docsrs

    lib_dir=$root_dir/pycmbs
    core_dir=$lib_dir/core
    bench_dir=$lib_dir/benchmarking
    conf_dir=$root_dir/configuration


    mkdir ${core_dir} -p
    mkdir ${bench_dir} -p
    mkdir ${conf_dir} -p

    touch $root_dir/pycmbs-benchmarking.py
    touch $lib_dir/__init__.py $bench_dir/__init__.py
    touch $core_dir/data.py $core_dir/grid.py $core_dir/__init__.py $core_dir/test_data.py $core_dir/test_grid.py
    touch $bench_dir/models.py $bench_dir/analysis.py $bench_dir/test_analysis.py $bench_dir/test_models.py
    touch $conf_dir/parameter.ini $conf_dir/models.json $conf_dir/analysis.json

    echo "Side-by-side layout"
    echo "###################"
    echo ""

    show_tree_output $root_dir
    rm -rf $root_dir
}

function nested_layout {

    echo ""
    echo "------"
    echo ""

    root_dir=/tmp/mock/pyCMBS-vX.YY
    mkdir $root_dir -p
    mkdir $root_dir/docs $root_dir/docsrs
    touch $root_dir/MANIFEST $root_dir/setup.py $root_dir/Makefile $root_dir/INSTALL \
            $root_dir/README $root_dir/LICENSE

    lib_dir=$root_dir/pycmbs
    core_dir=$lib_dir
    bench_dir=$lib_dir/benchmarking
    conf_dir=$root_dir/configuration


    mkdir ${core_dir} -p
    mkdir ${bench_dir} -p
    mkdir ${conf_dir} -p

    touch $root_dir/pycmbs-benchmarking.py
    touch $lib_dir/__init__.py $bench_dir/__init__.py
    touch $core_dir/data.py $core_dir/grid.py $core_dir/__init__.py $core_dir/test_data.py $core_dir/test_grid.py
    touch $bench_dir/models.py $bench_dir/analysis.py $bench_dir/test_analysis.py $bench_dir/test_models.py
    touch $conf_dir/parameter.ini $conf_dir/models.json $conf_dir/analysis.json

    echo "Nested layout"
    echo "#############"
    echo ""
    show_tree_output $root_dir
    rm -rf $root_dir
}


function single_lib_folder_layout {

    echo ""
    echo "------"
    echo ""

    root_dir=/tmp/mock/pyCMBS-vX.YY
    mkdir $root_dir -p
    mkdir $root_dir/docs $root_dir/docsrs
    touch $root_dir/MANIFEST $root_dir/setup.py $root_dir/Makefile $root_dir/INSTALL \
            $root_dir/README $root_dir/LICENSE

    lib_dir=$root_dir/pycmbs
    core_dir=$lib_dir
    bench_dir=$lib_dir
    conf_dir=$root_dir/configuration


    mkdir ${core_dir} -p
    mkdir ${bench_dir} -p
    mkdir ${conf_dir} -p

    touch $root_dir/pycmbs-benchmarking.py
    touch $lib_dir/__init__.py $bench_dir/__init__.py
    touch $core_dir/data.py $core_dir/grid.py $core_dir/__init__.py $core_dir/test_data.py $core_dir/test_grid.py
    touch $bench_dir/models.py $bench_dir/analysis.py $bench_dir/test_analysis.py $bench_dir/test_models.py
    touch $conf_dir/parameter.ini $conf_dir/models.json $conf_dir/analysis.json

    echo "Single lib folder (pycmbs) layout"
    echo "#############"
    echo ""
    show_tree_output $root_dir
    rm -rf $root_dir
}

# show layouts
side_by_side_layout
nested_layout
single_lib_folder_layout
