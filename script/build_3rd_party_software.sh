#/usr/bin/env bash
set -e
set -x


build_root=$HOME/tmp/lisem_external_build
install_prefix=$HOME/tmp/lisem_external


export PATH=/cygdrive/c/mingw64/bin:$PATH


function native_path()
{
    local pathname=$1
    local variable_name=$2

    native_pathname=`cygpath -m $pathname`
    eval $variable_name="$native_pathname"
}


unset PCRTEAM_EXTERN
cmake=cmake
cmake_generator="Unix Makefiles"
cmake_make_program=mingw32-make
wget=wget
unzip=unzip
sed="sed -i.tmp"  # Make sure this is GNU sed!


function build_using_cmake()
{
    local source_directory=$1
    local binary_directory=$2
    local install_prefix=$3
    local build_type=$4

    native_path $source_directory native_source_directory
    native_path $install_prefix native_install_prefix

    cd $binary_directory
    $cmake \
        -DCMAKE_BUILD_TYPE=$build_type \
        -G"$cmake_generator" \
        -DCMAKE_INSTALL_PREFIX="$native_install_prefix" \
        -DCMAKE_MAKE_PROGRAM="$cmake_make_program" \
        $native_source_directory
    $cmake --build . --config $build_type
    $cmake --build . --config $build_type --target install
}


function build_pcraster_raster_format()
{
    cd $build_root

    if [ ! -d pcraster_raster_format_sources ]; then
        rm -fr pcraster_raster_format_sources
        git clone git://git.code.sf.net/p/pcraster/rasterformat \
            pcraster_raster_format_sources
    fi

    rm -fr pcraster_raster_format_objects
    mkdir pcraster_raster_format_objects

    local raster_format_install_prefix=$install_prefix/pcraster_raster_format
    rm -fr $raster_format_install_prefix

    function build()
    {
        local build_type=$1
        build_using_cmake \
            $build_root/pcraster_raster_format_sources \
            $build_root/pcraster_raster_format_objects \
            $raster_format_install_prefix \
            $build_type
    }

    build Debug
    build Release
}


qt_version=4.8.4
read qt_main_version qt_sub_version qt_patch_version \
    < <(IFS=.; echo $qt_version)
qt_base=qt-everywhere-opensource-src-$qt_version
qt_make_spec=win32-g++  # Yes, 32.
qt_install_prefix=$install_prefix/qt-$qt_version


function build_qt()
{
    cd $build_root

    if [ ! -e $qt_base.zip ]; then
        $wget http://download.qt-project.org/archive/qt/$qt_main_version.$qt_sub_version/$qt_version/$qt_base.zip
    fi

    rm -fr $qt_base
    $unzip -n $qt_base.zip

    rm -fr $qt_install_prefix
    mv $qt_base $qt_install_prefix
    cd $qt_install_prefix

    # https://bugreports.qt-project.org/browse/QTBUG-34856
    # https://codereview.qt-project.org/#change,71532
    $sed "50s/ && !defined(Q_CC_GNU)//" src/corelib/tools/qsimd.cpp

    configure_arguments="
        -opensource
        -debug-and-release
        -no-qt3support
        -no-xmlpatterns
        -no-multimedia
        -no-phonon
        -no-phonon-backend
        -no-webkit
        -no-script
        -no-scripttools
        -no-declarative
        -nomake examples
        -nomake demos
    "
    accept_license_answer="y"
    ./configure.exe $configure_arguments <<ACCEPT_LICENSE
    $accept_license_answer
ACCEPT_LICENSE

    $cmake_make_program
}


build_qwt() {
    local qwt_version=6.0.1
    local qwt_base=qwt-$qwt_version
    local qwt_install_prefix=$install_prefix/$qwt_base

    cd $build_root

    if [ ! -e $qwt_base.tar.bz2 ]; then
        $wget http://downloads.sourceforge.net/project/qwt/qwt/$qwt_version/$qwt_base.tar.bz2
    fi

    rm -fr $qwt_base
    tar jxf $qwt_base.tar.bz2
    cd $qwt_base

    # Our install prefix.
    $sed "28s!^!QWT_INSTALL_PREFIX = $qwt_install_prefix!" qwtconfig.pri

    # We don't have svg support.
    $sed "s!QWT_CONFIG     += QwtSvg!# QWT_CONFIG     += QwtSvg!" qwtconfig.pri

    # # Uncomment an error message that is not relevant for us, but will stop
    # # the build.
    # $sed "57,69s/^/# /" qwtbuild.pri

    # http://ehc.ac/p/qwt/mailman/message/29258218/
    $sed "3625s/unsigned long/unsigned long long/" \
        textengines/mathml/qwt_mml_document.cpp

    $qt_install_prefix/bin/qmake -spec $qt_make_spec qwt.pro

    $cmake_make_program install
}


build_pcraster_raster_format
build_qt
build_qwt
