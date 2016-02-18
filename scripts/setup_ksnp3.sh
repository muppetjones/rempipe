#!/bin/bash
# Perform a local install of a variety of useful (e.g., necessary) softwares
#
# Ever had to use a server and found yourself frustrated due to lack of
# root access or simply not being able/allowed to modify the setup?
#   Never fear!
#
# This script does the following:
#   * Modifies *.vimrc to include line numbers and set tabs to four spaces
#   * Install gcc-4.5.4 (and other required)
#   * Install python-3.4 (also: venv called 'py3', activate replacement called 'start')
#   * Install bowtie2
#   * Install R-3.2.1
#   * Adds path/local/(lib|bin) to path via ~/.bash_profile (or ~/.profile)
#
# Also installs:
#   * kSNP3.0
#   * jellyfish-2.2.3
#   * MUMmer 3.23
#   * phylip-3.696
#
# To do:
#   * velvet, hisat
#   * test using CDPATH
#       export CDPATH=${root_prefix}:${local_prefix}
#   * finish moving download location to 'local/build' ($build_prefix)
#       -- can't yet...still testing on local machine, which has half in 'dev'.
#
# Sets up the following directory structure:
#   given/path/
#       dev/
#       local/
#           bin/
#           lib/
#       venv/

# exit on error!
set -e
set +x

#------------------------------------------------------------------------------
# Functions

pause(){
 read -n1 -rsp $'Press any key to continue or Ctrl+C to exit...\n'
}

speak() { echo -e "${NC}$*${NC}"; }
yell() { echo -e "\n$0: $*" >&2; }
die() { yell "ERROR: $*"; exit 111; }
try() { "$@" || die "cannot $*"; }

announce() {
    speak "${NC}------------------------------------------------------"
    speak $*
    pause
}

accept() {
    read -n1 -rsp $'Press \'y\' to continue or any key to exit: \n' response
    if [[ $response != 'y' ]]; then
        die "Cancelled by user"
    fi
}

check_if_line_exists()
{
    # grep wont care if one or both files dont exist.
    grep -qsFx "$LINE_TO_ADD" $FILE_ADD_TO
}

add_line_to_file()
{
    # profile=~/.profile
    # [ -w "$profile" ] || profile=~/.bash_profile
    echo "$LINE_TO_ADD" >> $FILE_ADD_TO
}

#------------------------------------------------------------------------------
# Cleanup!!
function finish {
    echo "Cleaning up..."
    for i in ${CLEANUP_FILES[@]}
    do
       echo "Removing file ${i}"
       rm -f "${i}"
    done
    for i in ${CLEANUP_DIRS[@]}
    do
       echo "Removing dir ${i}"
       rm -rf "${i}"
    done
}
trap finish EXIT
CLEANUP_FILES=()
CLEANUP_DIRS=()

function DO_NOT_CLEAN {
    keep="$1"
    CLEANUP_DIRS=( "${CLEANUP_DIRS[@]/%$keep}" )
    CLEANUP_FILES=( "${CLEANUP_FILES[@]/%$keep}" )
}

#------------------------------------------------------------------------------
# Colors

BLACK='\033[0;30m'
DARK_GRAY='\033[1;30m'
RED='\033[0;31m'
LIGHT_RED='\033[1;31m'
GREEN='\033[0;32m'
LIGHT_GREEN='\033[1;32m'
BROWN_ORANGE='\033[0;33m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
LIGHT_BLUE='\033[1;34m'
PURPLE='\033[0;35m'
LIGHT_PURPLE='\033[1;35m'
CYAN='\033[0;36m'
LIGHT_CYAN='\033[1;36m'
LIGHT_GRAY='\033[0;37m'
WHITE='\033[1;37m'

NC='\033[0m' # No Color
CY=${CYAN}  # highlight
GR=${GREEN}  # files / directories
PU=${PURPLE}  # executables
BO=${BROWN_ORANGE}  # misc

#------------------------------------------------------------------------------
# Setup

# set profile file
if [[ ! -s "$HOME/.bash_profile" && -s "$HOME/.profile" ]] ; then
  profile_file="$HOME/.profile"
else
  profile_file="$HOME/.bash_profile"
fi

# get the main install directory
root_prefix="/usr/opt/ "
if [ -z $1 ]; then
    if [ -d "$(pwd)/dev" ]; then
        root_prefix="$(pwd)"
    else  # ask the user
        speak "${NC}------------------------------------------------------"
        speak "Type the root directory and hit [ENTER]:${BROWN_ORANGE}"
        read root_prefix
        root_prefix="${root_prefix%%/}"  # remove trailing slash
    fi
fi

local_prefix="${root_prefix}/local"
share_prefix="${root_prefix}/share"
venv_prefix="${root_prefix}/venv"

dev_prefix="${local_prefix}/dev"
bin_prefix="${local_prefix}/bin"
lib_prefix="${local_prefix}/lib"
build_prefix="${local_prefix}/build"

# Alert the user
nl="\n  ${GR}"
speak "------------------------------------------------------"
speak "Software will be downloaded and compiled in:${nl}${build_prefix}"
speak "Virtual environments will be created in:${nl}${venv_prefix}"
speak "Binaries will be placed in:${nl}${bin_prefix}"
speak "Libraries will be place in:${nl}${lib_prefix}"
speak "*  NOTE: Some software downloads to:${nl}${dev_prefix}"
speak "** PATH will be updated via ${profile_file}"
speak "------------------------------------------------------"
accept

mkdir -p "${dev_prefix}"
mkdir -p "${venv_prefix}"
mkdir -p "${bin_prefix}"
mkdir -p "${lib_prefix}"
mkdir -p "${build_prefix}"
mkdir -p "${share_prefix}"


#------------------------------------------------------------------------------
# DO THESE FIRST!!

cd "${dev_prefix}"

########################
# Update PATH in profile file
if ! grep -q "PATH=${bin_prefix}:${lib_prefix}" "${profile_file}"; then
    announce "Adding PATH mods to ${GR}${profile_file}${NC} and sourcing..."
    FILE_ADD_TO="${profile_file}"
    LINE_TO_ADD="export PATH=${bin_prefix}:${lib_prefix}:\${PATH}"
    check_if_line_exists || add_line_to_file
fi

########################
# Update certs
if [ ! -d "${local_prefix}/certs" ] || ! grep -q "cacert.pem" ~/.gitconfig ;
then
    announce "Downloading updated certs to ${GR}${local_prefix}/certs${NC}"
    certsdir="${local_prefix}/certs"
    mkdir -p "${certsdir}"
    curl http://curl.haxx.se/ca/cacert.pem -o "${certsdir}/cacert.pem"

    gitconfig=~/.gitconfig
    cat >>"${gitconfig}" <<EOL
[http]
sslCAinfo = ${certsdir}/cacert.pem
EOL

fi

#------------------------------------------------------------------------------
# HELPERS

function get_downloaded_filename {
    local link=$1
    local f

    # OPTION 1: Sourceforge latest download
    if [[ "$link" == *sourceforge*latest*download* ]]; then
        f=$(LANG=C wget "$link" 2>&1 | sed -n "s/.*- \`\(.*\)' saved.*/\1/p")

    # OPTION 2: Anyplace else
    else
        f="$(basename "${link%%/download}")"
    fi

    echo "$f"
}

function get_root_compressed_directory {
    zip_file="$1"
    zip_cmd="$2"

    bdir="$(${zip_cmd} ${zip_file} | sed -e 's@^\./@@' | sed -e 's@/.*@@' | uniq)"
    echo "${bdir}"
}

function download_and_extract {
    # Download file(s) from a link and extract to desired directory

    local link=$1
    local __bdir
    cd "${build_prefix}";  # ensure we're in the build directory!

    if [[ "${link}" == *.git ]]; then
        __bdir="$(pwd)/$(basename "${link%%.git}")"
        if [ ! -d "${__bdir}" ]; then
            speak "Cloning from ${GR}${link}${NC}..."
            git clone "$link"
        fi
    else
        local zipf="$(pwd)/$(get_downloaded_filename "$link")"
        if [ ! -f "${zipf}" ]; then
            speak "Downloading '${GR}${zipf}${NC}'..."
            wget "${link}"  --no-check-certificate
        fi

        # determine compression type
        ext=${zipf##*\.}
        case "$ext" in
            gz|tgz)
                decomp="tar"
                decomp_opt="-zxf"
                decomp_list="tar -tzf"
                ;;
            bz2|tbz2)
                decomp="tar"
                decomp_opt="-jxf"
                decomp_list="tar -tjf"
                ;;
            tar)
                decomp="tar"
                decomp_opt="-xf"
                decomp_list="tar -tf"
                ;;
            zip)
                decomp="unzip"
                decomp_opt=""
                decomp_list="unzip -Z -1"
                ;;
            *)
                die "Unknown compression type"
                ;;
        esac

        # get the root directory
        # -- IFF there is more than one directory, use the file name
        __bdir="$(pwd)/$(get_root_compressed_directory "${zipf}" "${decomp_list}")"
        if [[ "$(wc -l <<< "$__bdir")" > 1 ]]; then
            # remove extension -- covers .tar*, .gz, .zip
            __bdir="${zipf%%.[tgz]*}"
            if [[ "$decomp" == tar ]]; then
                decomp="${decomp} -C ${__bdir}"
            else
                decomp="${decomp} -d ${__bdir}"
            fi
        fi

        # decompress
        speak "Extracting '${GR}${zipf}${NC}'..."
        ${decomp} ${decomp_opt} "${zipf}"

        # NOTE: Don't add to the cleanup list until we get to the end.
        #       Will prevent us from having to redownload if an error occurs
        #       while we're unzipping
        CLEANUP_FILES+=("${zipf}")
    fi

    # echo "${__bdir}";  # NOTE: it's absolute path
    # return the value!
    eval "$2='$__bdir'"
}

function download_and_link {
    # This function:
    #   1. downloads a compressed file from a given link into build dir
    #   2. decompresses the file and moves into the directory
    #   3. runs configuration and make
    #   4. links desired executables into bin dir
    #
    # Arguments:
    #   name (string): The name of the program (just for announcing)
    #   link (URL): URL of file download location.
    #   --config='<options'>: Options to pass to `./configure`
    #   --altdir=<str>: Provide if software extracts to a different directory.
    #   --link=exe,bin/exe: Comma-separated list of binaries to link to bin.
    #
    # NOTE: If no `--config` option is given, ALL leftover options will be
    #       Passed to `./configure`
    # NOTE: The `--config` option MUST be
    #       a) surrounded by double quotes
    #       b) have it's options surrounded by single quotes
    #       *) it's easier just to not use it
    #       Example)
    #           "--config='--enable-python-binding'"

    local name=$1
    local link=$2
    local opt="$*"
    local bdir=''

    cd "${build_prefix}";  # ensure we're in the build directory!
    announce "Building ${CY}${name}${NC} at '${GR}$(pwd)${NC}'"

    #### OPTIONAL parameters
    # build directory
    if [[ "${opt}" == *--renamedir* ]]; then
        local renamedir="${build_prefix}/${opt##*--renamedir=}"
        renamedir="${bdir%%[${IFS}]*}"
        opt="${opt//--altdir=*[${IFS}]/}"
    fi

    # link list
    # -- expects comma separated list
    if [[ "${opt}" == *--link* ]]; then
        local exec_link="${opt##*--link=}"
        exec_link="${exec_link%%[${IFS}]*}"
        opt="${opt//--link=*[${IFS}]/}"
    fi

    # configure options
    # expects exactly `--config='<options>'`
    local config
    if [[ "${opt}" == *--config* ]]; then
        config="${opt##*--config=\'}"
        config="${config%%\'*}"
    else
        config=
    fi

    #### START downloading, etc.
    # download (if we don't already have the expected directory)
    download_and_extract "$link" bdir
    echo "bdir: '${bdir}'"
    echo "pwd: '$(pwd)'"

    # link executables (if we were told to)
    if [ ! -z "${exec_link}" ]; then
        while IFS=',' read -ra EXEC_ARR; do
            for f in "${EXEC_ARR[@]}"; do
                speak "Linking ${PU}${f}${NC}..."
                f="${bdir}/${f}"
                if [ -f "${f}" ]; then
                    chmod +x "${f}"
                    ln -s "${f}" "${bin_prefix}"
                fi
            done
        done <<< "$exec_link"
    fi

    cd ..
}

function download_and_compile {
    # This function:
    #   1. downloads a compressed file from a given link into build dir
    #   2. decompresses the file and moves into the directory
    #   3. runs configuration and make
    #   4. links desired executables into bin dir
    #
    # Arguments:
    #   name (string): The name of the program (just for announcing)
    #   link (URL): URL of file download location.
    #   --config='<options'>: Options to pass to `./configure`
    #   --altdir=<str>: Provide if software extracts to a different directory.
    #   --link=exe,bin/exe: Comma-separated list of binaries to link to bin.
    #
    # NOTE: If no `--config` option is given, ALL leftover options will be
    #       Passed to `./configure`
    # NOTE: The `--config` option MUST be
    #       a) surrounded by double quotes
    #       b) have it's options surrounded by single quotes
    #       *) it's easier just to not use it
    #       Example)
    #           "--config='--enable-python-binding'"

    local name=$1
    local link=$2
    local opt="$*"
    local bdir=''

    cd "${build_prefix}";  # ensure we're in the build directory!
    announce "Building ${CY}${name}${NC} at '${GR}$(pwd)${NC}'"

    #### OPTIONAL parameters
    # build directory
    if [[ "${opt}" == *--renamedir* ]]; then
        local renamedir="${build_prefix}/${opt##*--renamedir=}"
        renamedir="${bdir%%[${IFS}]*}"
        opt="${opt//--altdir=*[${IFS}]/}"
    fi

    # link list
    # -- expects comma separated list
    if [[ "${opt}" == *--link* ]]; then
        local exec_link="${opt##*--link=}"
        exec_link="${exec_link%%[${IFS}]*}"
        opt="${opt//--link=*[${IFS}]/}"
    fi

    # configure options
    # expects exactly `--config='<options>'`
    local config
    if [[ "${opt}" == *--config* ]]; then
        config="${opt##*--config=\'}"
        config="${config%%\'*}"
    else
        config=
    fi

    #### START downloading, etc.
    # download (if we don't already have the expected directory)
    download_and_extract "$link" bdir
    echo "bdir: '${bdir}'"
    echo "pwd: '$(pwd)'"

    # rename downloaded directory if given
    if [ ! -z "${renamedir}" ]; then
        speak "Renaming '${GR}$(basename ${bdir})${NC}' to '${GR}${renamedir}${NC}"
        mv "${bdir}" "$(dirname "${bdir}")/${renamedir}"
        bdir="$(dirname "${bdir}")/${renamedir}"
    fi

    # configure and make (IFF we find a configure file)
    speak "Descending into ${GR}${bdir}${NC}"
    cd "${bdir}"
    if [ -f "${bdir}/configure" ]; then
        speak "Configuring and making..."
        ./configure --prefix="${local_prefix}" --exec_prefix="${local_prefix}" ${config}
        make
        make install
    else
        speak "Making..."
        make
    fi

    # link executables (if we were told to)
    if [ ! -z "${exec_link}" ]; then
        while IFS=',' read -ra EXEC_ARR; do
            for f in "${EXEC_ARR[@]}"; do
                speak "Linking ${PU}${f}${NC}..."
                f="${bdir}/${f}"
                if [ -f "${f}" ]; then
                    chmod +x "${f}"
                    ln -s "${f}" "${bin_prefix}"
                fi
            done
        done <<< "$exec_link"
    fi

    cd ..
}

## SIMPLE BUILD!
function build_from_tar_file() {
    cd "${build_prefix}"
    link="$1"
    tarfile="${build_prefix}/$(basename "$link")"
    dir="${tarfile%%.tar*}"

    CLEANUP_FILES+=("${tarfile}")
    wget "$link"
    tar zxvf "${tarfile}"
    cd "${dir}"

    # remove the first positional param ($link)
    shift
    ./configure --prefix="${local_prefix}" "$@"
    make
    make install
    cd "${build_prefix}"
}


#------------------------------------------------------------------------------
# Installs

########################
# g++ (into local!!)

# if [ ! -f "${local_prefix}/gcc-4.5.4/bin/gcc" ]; then
#     cd "${local_prefix}"
#
#     gmp="gmp-4.3.2"
#     mpfr="mpfr-2.4.2"
#     mpc="mpc-0.8.1"
#     gcc="gcc-4.5.4"
#     announce "Downloading and installing ${CY}${gcc}${NC}"
#     yell "-- thanks to solar via http://openwall.info/wiki/internal/gcc-local-build"
#
#     gmpdir="${local_prefix}/${gmp}"
#     mpfrdir="${local_prefix}/${mpfr}"
#     mpcdir="${local_prefix}/${mpc}"
#     gccdir="${local_prefix}/${gcc}"
#
#     mkdir -p "${gmpdir}"
#     mkdir -p "${mpfrdir}"
#     mkdir -p "${mpcdir}"
#     mkdir -p "${gccdir}"
#
#     speak "--- Downloading required tarballs into ${build_prefix}"
#     cd "${build_prefix}"
#     infrastlink="ftp://gcc.gnu.org/pub/gcc/infrastructure"
#     releaselink="ftp://gcc.gnu.org/pub/gcc/releases/${gcc}"
#
#     if [ ! -d "${gmpdir}/include" ]; then
#         speak "--- ${CY}${gmp}${NC}"
#         cd "${build_prefix}"
#         wget "${infrastlink}/${gmp}.tar.bz2"
#         tar xjf ${gmp}.tar.bz2
#         cd ${gmp}
#         ./configure --prefix="${gmpdir}" --enable-cxx
#         nice -n 19 time make -j8
#         make install
#         make check
#         echo $?
#         CLEANUP_FILES+=("${build_prefix}/${gmp}.tar.bz2")
#     fi
#
#     if [ ! -d "${mpfrdir}/include" ]; then
#         speak "--- ${CY}${mpfr}${NC}"
#         cd "${build_prefix}"
#         wget "${infrastlink}/${mpfr}.tar.bz2"
#         tar xjf ${mpfr}.tar.bz2
#         cd ${mpfr}
#         ./configure --prefix="${mpfrdir}" --with-gmp="${gmpdir}"
#         nice -n 19 time make -j8
#         make install
#         CLEANUP_FILES+=("${build_prefix}/${mpfr}.tar.bz2")
#     fi
#
#     if [ ! -d "${mpcdir}/include" ]; then
#         speak "--- ${CY}${mpc}${NC}"
#         cd "${build_prefix}"
#         wget "${infrastlink}/${mpc}.tar.gz"
#         tar xzf ${mpc}.tar.gz
#         cd ${mpc}
#         LD_LIBRARY_PATH=${gmpdir}/lib:${mpfrdir}/lib ./configure --prefix=${mpcdir} --with-gmp=${gmpdir} --with-mpfr=${mpfrdir}
#         LD_LIBRARY_PATH=${gmpdir}/lib:${mpfrdir}/lib nice -n 19 time make -j8
#         make install
#         CLEANUP_FILES+=("${build_prefix}/${mpc}.tar.gz")
#     fi
#
#     if [ ! -d "${gccdir}/include" ]; then
#         speak "--- ${CY}${gcc}${NC}"
#         cd "${build_prefix}"
#         if [ ! -d "${build_prefix}/${gcc}" ]; then
#             wget "${releaselink}/${gcc}.tar.bz2"
#             tar xjf ${gcc}.tar.bz2
#         fi
#         cd "${build_prefix}/${gcc}"
#         LD_LIBRARY_PATH=${gmpdir}/lib:${mpfrdir}/lib:${mpcdir}/lib ./configure --prefix=${gccdir} --with-gmp=${gmpdir} --with-mpfr=${mpfrdir} --with-mpc=${mpcdir} --host=x86_64-redhat-linux --enable-languages=c,c++,fortran --disable-multilib --enable-threads=posix
#         LD_LIBRARY_PATH=${gmpdir}/lib:${mpfrdir}/lib:${mpcdir}/lib nice -n 19 time make -j8
#         make install
#         CLEANUP_FILES+=("${build_prefix}/${gcc}.tar.bz2")
#     fi
#
#     FILE_ADD_TO="${profile_file}"
#     LINE_TO_ADD="export LD_LIBRARY_PATH=${gmpdir}/lib:${mpfrdir}/lib:${mpcdir}/lib:${gccdir}/lib64:\${LD_LIBRARY_PATH}"
#     check_if_line_exists || add_line_to_file
#     LINE_TO_ADD="export PATH=${gccdir}/bin:\$PATH"
#     check_if_line_exists || add_line_to_file
#
# fi
#
#
# ########################
# # autotools -- autconf and automake
# # NOTE: current system has m4 vX.XX, which works fine
# if  [ ! -f "${bin_prefix}/autoconf" ] || [ ! -f "${bin_prefix}/automake" ];
# then
#     if [ ! -f "${bin_prefix}/pkg-config" ]; then
#         announce "Downloading and compiling ${CY}pkg-config${NC}"
#         link="http://pkgconfig.freedesktop.org/releases/pkg-config-0.28.tar.gz"
#         build_from_tar_file "$link"
#         hash pkg-config
#     fi
#
#     if [ ! -f "${bin_prefix}/m4" ]; then
#         announce "Downloading and compiling ${CY}m4${NC}"
#         link="http://ftp.gnu.org/gnu/m4/m4-latest.tar.gz"
#         link="http://ftp.gnu.org/gnu/m4/m4-1.4.17.tar.gz"
#
#         # NOTE: The change in CPPFLAGS is necessary for older systems
#         build_from_tar_file "$link" CPPFLAGS="-fgnu89-inline"
#         hash m4
#     fi
#
#     if [ ! -f "${bin_prefix}/libtool" ]; then
#         announce "Downloading and compiling ${CY}libtool${NC}"
#         link="http://ftp.gnu.org/gnu/libtool/libtool-latest.tar.gz"
#         link="http://ftp.gnu.org/gnu/libtool/libtool-2.4.6.tar.gz"
#         build_from_tar_file "$link"
#         hash libtool
#         hash libtoolize
#     fi
#
#     if [ ! -f "${bin_prefix}/autoconf" ]; then
#         announce "Downloading and compiling ${CY}autoconf${NC}"
#         link="http://ftp.gnu.org/gnu/autoconf/autoconf-latest.tar.gz"
#         link="http://ftp.gnu.org/gnu/autoconf/autoconf-2.69.tar.gz"
#         link="http://ftp.gnu.org/gnu/autoconf/autoconf-2.65.tar.gz"
#         build_from_tar_file "$link"
#         hash autoconf
#         hash autoreconf
#         hash autoheader
#         hash autom4te
#         hash autoscan
#         hash autoupdate
#     fi
#     if [ ! -f "${bin_prefix}/automake" ]; then
#         announce "Downloading and compiling ${CY}automake${NC}"
#         link="http://ftp.gnu.org/gnu/automake/automake-latest.tar.gz"
#         link="http://ftp.gnu.org/gnu/automake/automake-1.15.tar.gz"
#         build_from_tar_file "$link"
#         hash automake
#         hash aclocal
#     fi
# fi
#
#
# ########################
# # python3
# if [ ! -f "${bin_prefix}/python3" ]; then
#     announce "Setting up ${CY}Python-3.4.3${NC}..."
#     cd "${build_prefix}"
#
#     link="https://www.python.org/ftp/python/3.4.3/Python-3.4.3.tgz"
#     pytar="$(basename ${link})"
#     pydir="${build_prefix}/$(basename ${link%%.tgz})"
#
#     CLEANUP_DIRS+=("${pydir}")
#     CLEANUP_FILES+=("${pytar}")
#     CLEANUP_DIRS+=("${venv_prefix}/py3")
#     CLEANUP_FILES+=("${bin_prefix}/python3")
#     CLEANUP_FILES+=("${bin_prefix}/pip3")
#     CLEANUP_FILES+=("${bin_prefix}/pyvenv")
#
#     speak "Downloading and compiling ${GR}${pytar}${NC}..."
#     wget ${link} --no-check-certificate
#     tar zxvf Python-3.4.3.tgz
#     cd Python-3.4.3
#
#     ./configure --prefix="${local_prefix}" --exec_prefix="${local_prefix}"
#     make
#     make altinstall
#
#     # links
#     speak "Linking ${PU}python3${NC}, ${PU}pip3${NC}, ${PU}pyenv${NC}"
#     ln -s "${bin_prefix}/python3.4" "${bin_prefix}/python3"
#     ln -s "${bin_prefix}/pip3.4" "${bin_prefix}/pip3"
#     ln -s "${bin_prefix}/pyvenv-3.4" "${bin_prefix}/pyvenv"
#
#     # venv start script
#     speak "Creating virtual environment script ${PU}start${NC}"
#     scriptfile="${bin_prefix}/start"
#     cat >"${scriptfile}" <<EOL
# #!/bin/bash
# # A replacement for activate, based (shamelessly copied) on 'inve'
# # -- thanks to datagrok via https://gist.github.com/datagrok/2199506
# # Improvements (modifications?) for 'deactivate' and .bashrc from muppetjones
# # NOTE: If you add to the front of your path in your .bashrc
# #       or .bash_profile, it will overwrite these settings.
# #       Consider using a conditional setup for VIRTUAL_ENV.
# #           if [ -z "$VIRTUAL_ENV" ]; then
# #               export PATH=/mypath/bin:${PATH}
# #           else
# #               alias deactivate='exit'
# #               export PATH=${PATH}:/mypath/bin:${lib_path}
# #           fi
# # NOTE: To include prompt replacement, add the following to your .bashrc, etc:
# #           __venv_ps1() {
# #               if [ ! -z $VIRTUAL_ENV ]; then
# # 		            printf -- "(%s) " "$(basename ${VIRTUAL_ENV})"
# # 	            fi
# #           }
# #           PS1='$(declare -F __venv_ps1 &>/dev/null && __venv_ps1)$PS1'
#
# export VIRTUAL_ENV="${venv_prefix}/\$1"
# export PATH="$VIRTUAL_ENV/bin:$PATH"
# unset PYTHONHOME
# exec "$SHELL"
#
# # Just to keep with tradition, allow the use of 'deactivate'
# # NOTE: You must set "alias deactivate='exit'" in your .bashrc
# #       or .bash_profile, preferably conditionally defined
# #       within [ ! -z ${VIRTUAL_ENV} ] (see above example)
#
# unalias deactivate
# EOL
#
#     chmod +x "${scriptfile}"
#
#     speak "Creating virtual env ${BO}py3${NC} in ${GR}${venv_prefix}${NC}"
#     cd "${venv_prefix}"
#     pyvenv py3
#
#     DO_NOT_CLEAN "${dev_prefix}/Python-3.4.3"
#     DO_NOT_CLEAN "${venv_prefix}/py3"
#     DO_NOT_CLEAN "${dev_prefix}/Python-3.4.3.tgz"
#     DO_NOT_CLEAN "${bin_prefix}/python3"
#     DO_NOT_CLEAN "${bin_prefix}/pip3"
#     DO_NOT_CLEAN "${bin_prefix}/pyvenv"
#
#     cd "${dev_prefix}"
# fi
#
# ########################
# # bowtie2
# if [ ! -f "${bin_prefix}/bowtie2" ]; then
#     name="bowtie2"
#     link="git@github.com:BenLangmead/bowtie2.git"
#     bin_list="bowtie2,bowtie2-build,bowtie2-inspect"
#
#     download_and_compile $name $link --link=$bin_list
# fi
#
# ########################
# # R
#
# if [ ! -f "${bin_prefix}/R" ]; then
#     name="R"
#     link="https://cran.r-project.org/src/base/R-3/R-3.2.1.tar.gz"
#
#     # NOTE: R requires gcc to be compiled for fortran
#     download_and_compile $name $link
# fi
#
#
# ########################
# # fastqc
#
# if [ ! -f "${bin_prefix}/fastqc" ]; then
#     name="FastQC"
#     link="http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.3.zip"
#     altdir="FastQC"
#     bin_list="fastqc"
#
#     download_and_compile $name $link --link=$bin_list
# fi
#
#
# ########################
# # TagDust
#
# if [ ! -f "${bin_prefix}/tagdust" ]; then
#     name="TagDust"
#     link="http://sourceforge.net/projects/tagdust/files/tagdust-2.31.tar.gz/download"
#     bin_list="tagdust2"
#
#     download_and_compile $name $link --link=$bin_list
# fi
#
#
# ########################
# # Skewer
#
# if [ ! -f "${bin_prefix}/skewer" ]; then
#     name="Skewer"
#     link="http://sourceforge.net/projects/skewer/files/Binaries/skewer-0.1.127-linux-x86_64/download"
#     bin_list="skewer"
#
#     cd "${bin_prefix}"
#     wget "${link}"
#     ln -s "${bin_prefix}/skewer-0.1.127-linux-x86_64" "skewer"
#     chmod +x skewer*
# fi
#
# ########################
# # GNU scientific library
#
# if [ ! -f "${bin_prefix}/gsl-1.6"* ]; then
#     name="GSL"
#     link="http://gnu.mirror.vexxhost.com/gsl/gsl-1.16.tar.gz"
#
#     download_and_compile $name $link
#
#     FILE_ADD_TO="${profile_file}"
#     LINE_TO_ADD="export LD_LIBRARY_PATH=${lib_prefix}:${lib_prefix}64:\${LD_LIBRARY_PATH}"
#     check_if_line_exists || add_line_to_file
#
#     FILE_ADD_TO="${profile_file}"
#     LINE_TO_ADD="export CPATH=${local_prefix}/include:\${CPATH}"
#     check_if_line_exists || add_line_to_file
# fi
#
# ########################
# # EA-utils
#
# if [ ! -f "${bin_prefix}/fastq-mcf" ]; then
#     name="EA-utils"
#     announce "Downloading and installing ${CY}${name}${NC}"
#
#     link="http://ea-utils.googlecode.com/svn/trunk/"
#     bin_list="skewer"
#
#     cd "${build_prefix}"
#     svn checkout http://ea-utils.googlecode.com/svn/trunk/ ea-utils-read-only
#     cd ea-utils-read-only/clipper
#     PREFIX="${local_prefix}" make install
#
#
# fi


########################
# kSNP

# if [ ! -f "${bin_prefix}/kSNP3" ]; then
#     announce "Downloading and installing ${CY}kSNP 3.0${NC}"
#     cd "${dev_prefix}"
#
#     link="http://sourceforge.net/projects/ksnp/files/kSNP3.0_Linux_package.zip/download"
#     ksnptar="$(basename $(dirname "${link}"))"
#     ksnpdir="${dev_prefix}/${ksnptar%%.zip}"
#     ksnpdir_rename="kSNP3_Linux"
#
#     CLEANUP_FILES+=("${dev_prefix}/${ksnptar}")
#
#     if [ ! -d "${ksnpdir_rename}" ]; then
#         wget "${link}" --no-check-certificate
#         unzip "${ksnptar}"
#         mv "${ksnpdir}" "kSNP3_Linux"
#         perl -pi -e "s|(set kSNP=)/usr/local/kSNP3|\1$(pwd)/kSNP3|" "${ksnpdir}/kSNP3/kSNP3"
#     fi
#     ksnpdir="${dev_prefix}/kSNP3_Linux"
#     cd "${ksnpdir}"
#
#
#     cat >"${ksnpdir}/run_me.sh" <<EOL
# #!/bin/bash
# # Run __prog__ in whatever folder you darn well please.
# PATH=${ksnpdir}/kSNP3:$PATH
# "${ksnpdir}/kSNP3/__prog__" "\$@"
# EOL
#
#     function make_ksnp_command() {
#         prog_name=$1
#         script_name="${ksnpdir}/kSNP3/run_${prog_name}.sh"
#         cp "${ksnpdir}/run_me.sh" "${script_name}"
#         perl -pi -e "s/__prog__/${prog_name}/g" "${script_name}"
#         ln -s "${script_name}" "${bin_prefix}/${prog_name}"
#         chmod +x "${script_name}"
#     }
#
#     make_ksnp_command Kchooser
#     make_ksnp_command kSNP3
#     make_ksnp_command MakeFasta
#
# fi


########################
# jellyfish -- BREWABLE

# if [ ! -d "${dev_prefix}/jellyfish-2.2.3" ]; then
#     announce "Cloning and installing ${CY}Jellyfish${NC}"
#     cd "${build_prefix}"
#
#     # github install
#     # -- requires g++ >= 4.4
#     # link="git@github.com:gmarcais/Jellyfish.git"
#     # gdir="${dev_prefix}/$(basename "${gitssh%%.git}")"
#     # CLEANUP_DIRS+=("${gdir}")
#     # git clone "${link}"
#     # cd "${gdir}"
#     # autoreconf -i
#     # ./configure --prefix="${local_prefix}" --enable-python-binding="${lib_prefix}/python3.4"
#     # make
#     # make install
#
#     link="https://github.com/gmarcais/Jellyfish/releases/download/v2.2.3/jellyfish-2.2.3.tar.gz"
#     tarfile="$(basename "$link")"
#     dir="$(pwd)/${tarfile%%.tar.gz}"
#     CLEANUP_DIRS+=("${dir}")
#     CLEANUP_FILES+=("${build_prefix}/${tarfile}")
#
#     curl -O -L -k "${link}"  # redirects to amazon s3...very frustrating
#     tar zxvf "${tarfile}"
#
#     cd "${dir}"
#     ./configure --prefix="${local_prefix}" --enable-python-binding="${lib_prefix}/python3.4"
#     make
#     make install
#
#     DO_NOT_CLEAN "${dir}"
#     cd "${dev_prefix}"
# fi

########################
# Parsimonator

if [ ! -f "${bin_prefix}/parsimonator" ]; then
    announce "Cloning and installing ${CY}Parsimonator${NC}"
    cd "${build_prefix}"
    link="git@github.com:stamatak/Parsimonator-1.0.2.git"
    gdir="$(pwd)/$(basename "${link%%.git}")"
    git clone "${link}"

    cd "${gdir}"
    make -f Makefile.gcc

    ln -s "${gdir}/parsimonator" "${bin_prefix}/"
fi

########################
# Mummer -- BREWABLE
#
# if [ ! -f "${bin_prefix}/mummer" ]; then
#     announce "Downloading and installing ${CY}Mummer${NC}"
#     cd "${build_prefix}"
#     link="http://sourceforge.net/projects/mummer/files/latest/download"
#     file=$(LANG=C wget "${link}" 2>&1 | sed -n "s/.*- \`\(.*\)' saved.*/\1/p")
#     dir="$(pwd)/$(basename "${file%%.tar.gz*}")"
#
#     CLEANUP_FILES+=("${build_prefix}/${file}")
#
#     tar zxvf "${file}"
#     cd "${dir}"
#
#     make check
#     make install
#
#     for f in *; do
#         if [-f "${f}" ] && [ -x "${f}" ]; then
#             ln -s "$(pwd)/${f}" "${bin_prefix}/"
#         fi
#     done
# fi

########################
# Primux -- already in kSNP (and a pain to get working)
# if [ ! -f "${bin_prefix}/primux" ]; then
#     announce "Downloading and installing ${CY}Primux${NC}"
#     # NOTE: The tar creates a 'primux' directory (no version name)
#
#     cd "${build_prefix}"
#     set -x
#     link="http://sourceforge.net/projects/primux/files/primux_20_july_2014.tar.gz/download"
#     file="$(basename $(dirname "${link}"))"
#     dir="$(pwd)/${filename%%.tar*}"
#     wget -O "${file}" "${link}"
#
#     CLEANUP_FILES+=("${file}")
#
#     tar zxvf "${file}"
#     cd "${dir}"
#
#     make check
#     make install
# fi


########################
# Phylip -- BREWABLE
# -- kSNP uses consense, but it is finicky
# if [ ! -f "${build_prefix}/phylip-3.696/exe/consense" ]; then
#     announce "Downloading and installing ${CY}Phylip${NC}"
#     cd "${build_prefix}"
#     link="http://evolution.gs.washington.edu/phylip/download/phylip-3.696.tar.gz"
#     file="$(basename "${link}")"
#     dir="$(pwd)/${file%%.tar*}"
#
#     CLEANUP_FILES+=("${build_prefix}/${file}")
#
#     if [ ! -d "$dir/src" ]; then
#         wget "$link"
#         tar zxvf "$file"
#         cd "$dir/src"
#     fi
#
#     if [ ! -f "$dir/exe/consense" ]; then
#         make -f Makefile.unx install
#     fi
#
#     # fix the kSNP situation
#     ksnplink="$(readlink -f "${bin_prefix}/kSNP3")"
#     ksnpdir="$(dirname "${ksnplink}")"
#     if [ -d "${ksnpdir}" ]; then
#         mkdir -p "${ksnpdir}/backup"
#         mv "${ksnpdir}/consense" "${ksnpdir}/backup/"
#         ln -s "${dir}/exe/consense" "${ksnpdir}"
#     fi
# fi

########################
# RAST (and SEED)
# -- annotate prokaryotic genomes
# -- does NOT contain parent directory!!
# if [[ ! "${PATH}" !=  *"/sas/bin"* ]]; then
#     announce "Downloading and installing ${CY}RAST${NC}"
#     cd "${build_prefix}"
#     link="http://blog.theseed.org/downloads/sas-33a23.tgz"
#     link="http://blog.theseed.org/downloads/sas.tgz"
#     file="$(basename "${link}")"
#     dir="$(pwd)/${file%%.tgz}"
#
#     CLEANUP_FILES+=("${dir}/${file}")
#
#     if [ ! -d "$dir/modules" ]; then
#         mkdir -p "${dir}"
#         cd "${dir}"
#         wget "$link"
#         tar zxvf "$file"
#     fi
#
#     cd "${dir}/modules"
#     ./BUILD_MODULES
#
#     FILE_ADD_TO="${profile_file}"
#     LINE_TO_ADD="export PERL5LIB=\$PERL5LIB:${dir}/lib:${dir}/modules/lib"
#     check_if_line_exists || add_line_to_file
#     LINE_TO_ADD="export PATH=\$PATH:${dir}/bin"
#     check_if_line_exists || add_line_to_file
# fi

#######################
# Pstat
# if [ ! -f "${bin_prefix}/pstat" ]; then
#     cd "${build_prefix}"
#     link="git@github.com:muppetjones/qstat-pretty.git"
#     gdir="$(pwd)/$(basename "${link%%.git}")"
#     git clone "${link}"
#
#     cd "${gdir}"
#     python setup.py build
#     ln -s "${gdir}/pstat" "${bin_prefix}/"
# fi

# source profile
source "${profile_file}"
