#!/bin/bash
show_help() {
    echo "Usage:"
    echo "$(basename $0) host-notebook --port 8888 --listen 127.0.0.1"
    echo "$(basename $0) run-tests"
}
if [ $# -eq 0 ]; then
    show_help
    exit 1
fi
set -euxo pipefail

HOST_NOTEBOOK=0
RUN_TESTS=0
PORT=8000
while [ $# -gt 0 ]; do
    case "$1" in
        host-notebook)
            HOST_NOTEBOOK=1
            shift
            ;;
        --port)
            shift
            PORT=$1
            shift
            ;;
        --listen)
            shift
            LISTEN=$1
            shift
            ;;
        run-tests)
            RUN_TESTS=1
            shift
            ;;
        --help|-h)
            show_help
            exit 0
            ;;
        --)
            break;
            ;;
        *)
            show_help
            exit 1
            ;;
    esac
done

if ! which podrun >/dev/null; then
    SOME_TEMP_DIR=$(mktemp -d)
    trap 'rm -rf -- "$SOME_TEMP_DIR"' EXIT
    ( cd "$SOME_TEMP_DIR"; curl -LOs https://raw.githubusercontent.com/bjodah/dotfiles/master/per-leaf-dir/bin/podrun; chmod +x podrun )
    export PATH="$SOME_TEMP_DIR:$PATH"
fi

SCRIPTS_DIR=$(dirname $(realpath "$BASH_SOURCE"))
REPO_DIR=$(realpath "$SCRIPTS_DIR/..")
ENV_DIR="$REPO_DIR/.env"
PKG_NAME=${PKG_NAME:-$(basename $REPO_DIR)}

mkdir -p $ENV_DIR

cat <<EOF>$ENV_DIR/setup.sh
if [ ! -d .env/ ]; then
   >&2 echo "No .env directory?"
   exit 1
fi
if [ ! -d .env/venv ]; then
   python3 -m venv .env/venv   
fi
source .env/venv/bin/activate
if ! python -c "import pycvodes" 2>&1 >/dev/null; then
    python -m pip install --upgrade-strategy=eager --upgrade pip && \
    python -m pip install --upgrade-strategy=eager numpy 'cython>=3.0.10' setuptools && \
    env \
        PYCVODES_NO_LAPACK=1 \
        PYCVODES_NO_KLU=1 \
        LDFLAGS='-Wl,--disable-new-dtags -Wl,-rpath,/usr/local/lib -L/usr/local/lib' \
        python -m pip install --cache-dir .env/pypi-cache -e .[all]
fi
if ! python -c "import pyodesys" 2>&1 >/dev/null; then
   python -m pip install --cache-dir .env/pypi-cache -e .[all]
   #jupyter-nbextension enable --user --py widgetsnbextension
fi
export MPLBACKEND=Agg
EOF

cat <<EOF>$ENV_DIR/run-tests.sh
#!/bin/bash
set -e
source .env/setup.sh
pytest -sv -ra --pyargs $PKG_NAME
EOF

cat <<EOF>$ENV_DIR/host-notebook.sh
#!/bin/bash
set -e
source .env/setup.sh
jupyter notebook --no-browser --port $PORT --ip=* index.ipynb
EOF


if [ $RUN_TESTS = 1 ]; then
    podrun --cont-img-dir $SCRIPTS_DIR/environment \
           --name "${PKG_NAME}_run_tests" \
           -- bash $ENV_DIR/run-tests.sh
fi
if [ $HOST_NOTEBOOK = 1 ]; then
    podrun --cont-img-dir $SCRIPTS_DIR/environment \
           --name "${PKG_NAME}_host_notebook_${PORT}" \
           -p $LISTEN:$PORT:$PORT \
           -- bash $ENV_DIR/host-notebook.sh
fi
