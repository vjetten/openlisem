#/usr/bin/env bash
set -e
set -x


unset PCRTEAM_EXTERN


# TODO Add arguments:
# - build_root
# - install_prefix
build_root=$HOME/tmp/lisem_external_build
install_prefix=$HOME/tmp/lisem_external

mkdir -p $build_root $install_prefix


if [ `uname -o 2>/dev/null` ]; then
    os=`uname -o`
else
    os=`uname`
fi


if [ $os == "Cygwin" ]; then
    # Path to ml64.exe, required by Boost.Context.
    vs_2008_root=`cygpath "$VS90COMNTOOLS"`
    amd64_root="$vs_2008_root/../../VC/BIN/amd64"

    # Path to compiler.
    mingw_root=/cygdrive/c/mingw64/bin

    export PATH="$mingw_root:$amd64_root:$PATH"

    make=mingw32-make
else
    make=make
fi


cmake=cmake
cmake_generator="Unix Makefiles"
find=/usr/bin/find
cmake_make_program=$make
wget=wget
unzip=unzip
sed="sed -i.tmp"  # Make sure this is GNU sed!


raster_format_version=head
raster_format_install_prefix=$install_prefix/pcraster_raster_format-$raster_format_version

qt_version=4.8.4
read qt_main_version qt_sub_version qt_patch_version \
    < <(IFS=.; echo $qt_version)
qt_base=qt-everywhere-opensource-src-$qt_version
qt_make_spec=win32-g++  # Yes, 32.
qt_install_prefix=$install_prefix/qt-$qt_version

qwt_version=6.0.1
qwt_base=qwt-$qwt_version
qwt_install_prefix=$install_prefix/$qwt_base

boost_version=1.55.0
boost_base=boost_${boost_version//./_}
boost_install_prefix=$install_prefix/boost-${boost_version}

fern_version=head
fern_install_prefix=$install_prefix/fern-$fern_version


function native_path()
{
    local pathname=$1
    local variable_name=$2
    local native_pathname=$pathname

    if [ $os == "Cygwin" ]; then
        native_pathname=`cygpath -m $native_pathname`
    fi

    eval $variable_name="$native_pathname"
}


function rebuild_using_cmake()
{
    local source_directory=$1
    local binary_directory=$2
    local install_prefix=$3
    local build_type=$4
    local options=$5

    native_path $source_directory native_source_directory
    native_path $install_prefix native_install_prefix

    mkdir -p $binary_directory
    cd $binary_directory
    rm -fr *
    $cmake \
        -DCMAKE_BUILD_TYPE=$build_type \
        -G"$cmake_generator" \
        -DCMAKE_MAKE_PROGRAM=$cmake_make_program \
        -DCMAKE_INSTALL_PREFIX="$native_install_prefix" \
        $options \
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

    rm -fr $raster_format_install_prefix

    function build()
    {
        local build_type=$1
        rebuild_using_cmake \
            $build_root/pcraster_raster_format_sources \
            $build_root/pcraster_raster_format_objects \
            $raster_format_install_prefix \
            $build_type
    }

    if [ $os == "Cygwin" ]; then
        build Debug
    fi

    build Release
}


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
        -platform $qt_make_spec
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

    $make
    $make clean
    rm -fr demos examples qmake tmp
    $find . -name "*.o" -exec rm -f {} \;

    # This will make the qt installation moveable.
    echo -e "[Paths]\nPrefix=.." > bin/qt.conf
}


function build_qwt() {
    cd $build_root

    if [ ! -e $qwt_base.tar.bz2 ]; then
        $wget http://downloads.sourceforge.net/project/qwt/qwt/$qwt_version/$qwt_base.tar.bz2
    fi

    rm -fr $qwt_base
    tar jxf $qwt_base.tar.bz2
    cd $qwt_base

    # Our install prefix.
    $sed "28s!^!QWT_INSTALL_PREFIX = $qwt_install_prefix!" qwtconfig.pri

    ### # We don't have svg support.
    ### $sed "s!QWT_CONFIG     += QwtSvg!# QWT_CONFIG     += QwtSvg!" qwtconfig.pri

    # http://ehc.ac/p/qwt/mailman/message/29258218/
    $sed "3625s/unsigned long/unsigned long long/" \
        textengines/mathml/qwt_mml_document.cpp

    $qt_install_prefix/bin/qmake -spec $qt_make_spec qwt.pro

    $make install
}


function build_boost() {
    local boost_toolset=mingw
    local boost_address_model=64
    local boost_variant="debug,release"
    # local boost_link="static,shared"
    local boost_link="shared"
    local boost_threading=multi

    cd $build_root

    if [ ! -e $boost_base.tar.bz2 ]; then
        $wget http://downloads.sourceforge.net/project/boost/boost/$boost_version/$boost_base.tar.bz2
    fi

    rm -fr $boost_base
    tar jxf $boost_base.tar.bz2

    cd $boost_base
    cd tools/build/v2
    ./bootstrap.sh --with-toolset=$boost_toolset
    cp engine/bin.ntx86_64/{b2,bjam}.exe ../../..
    cd ../../..
    rm -fr $boost_install_prefix
    # http://www.boost.org/boost-build2/doc/html/bbv2/overview/invocation.html
    ./b2 \
        -j 4 \
        --layout=tagged \
        --prefix=`cygpath -m $boost_install_prefix` \
        --without-mpi \
        --without-python \
        toolset=gcc \
        variant=$boost_variant \
        address-model=$boost_address_model \
        link=$boost_link \
        threading=$boost_threading \
        install
}


function build_fern()
{
    cd $build_root

    if [ ! -d fern_sources ]; then
        rm -fr fern_sources
        git clone ssh://root@agung/mnt/Depot/fern fern_sources
        cd fern_sources
        git checkout feature/operation/laplacian
        cd ..
    fi

    rm -fr $fern_install_prefix

    function build()
    {
        local build_type=$1
        local options="
            -DFERN_ALGORITHM_ONLY:BOOL=TRUE
            -DBOOST_ROOT=`cygpath -m $boost_install_prefix`
        "
        rebuild_using_cmake \
            $build_root/fern_sources \
            $build_root/fern_objects \
            $fern_install_prefix \
            $build_type \
            "$options"

        # TODO Test fails in release build.
        # PATH="$boost_install_prefix/lib:$PATH" $cmake \
        #     --build `cygpath -m $build_root/fern_objects` \
        #     --config $build_type \
        #     --target test
    }

    if [ $os == "Cygwin" ]; then
        build Debug
    fi

    build Release
}


function create_zip()
{
    install_prefix_basename=`basename $install_prefix`
    date=`date +%Y%m%d`
    zip_filename=${install_prefix_basename}-${date}.zip

    cd $install_prefix/..
    rm -f $zip_filename
    zip -r -q -9 $zip_filename $install_prefix_basename
    echo `pwd`/$zip_filename
}


build_pcraster_raster_format
build_qt
build_qwt
build_boost
build_fern
create_zip
