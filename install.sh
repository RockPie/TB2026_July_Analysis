#!/bin/sh

set -e

if ! command -v sudo > /dev/null; then
    echo "The sudo command is not installed. Please install it before running this script."
    exit 1
fi

ARCH="$(uname -m)"

if [ ! "$ARCH" = "x86_64" ] && [ ! "$ARCH" = "amd64" ] && [ ! "$ARCH" = "aarch64" ] && [ ! "$ARCH" = "arm64" ]; then
    printf "Your architecture: %s, is not supported." "$ARCH"
    exit 1
fi

. "/etc/os-release"

install_deb() {
    if [ ! "$ARCH" = "x86_64" ] && [ ! "$ARCH" = "amd64" ]; then
	    printf "Your architecture: %s, is not supported.\nYour architecture might, however, be supported by the Pip package. You may try the following instructions: https://docs.eduvpn.org/client/linux/installation.html#pip-installation\n" "$ARCH"
	    exit 1
    fi

    set -xeuo pipefail
    sudo apt-get update
    # Make sure dependencies are installed:
    # * apt-transport-https required to fetch packages over https
    # * curl + ca-certificates required to fetch the key
    # * gnupg required to verify the key
    sudo apt-get install --no-install-recommends apt-transport-https curl gnupg ca-certificates
    curl -sSf https://app.eduvpn.org/linux/v4/deb/app+linux@eduvpn.org.asc | gpg --dearmor | sudo tee /usr/share/keyrings/eduvpn-v4.gpg > /dev/null
    echo "deb [arch=$(dpkg --print-architecture) signed-by=/usr/share/keyrings/eduvpn-v4.gpg] https://app.eduvpn.org/linux/v4/deb/ $1 main" | sudo tee /etc/apt/sources.list.d/eduvpn-v4.list
    sudo apt-get update
    sudo apt-get install eduvpn-client
    exit 0
}

install_fedora() {
    set -x
    curl https://app.eduvpn.org/linux/v4/rpm/app+linux@eduvpn.org.asc | sudo tee /etc/pki/rpm-gpg/RPM-GPG-KEY-python-eduvpn-client_v4 > /dev/null
    cat << 'EOF' | sudo tee /etc/yum.repos.d/python-eduvpn-client_v4.repo
[python-eduvpn-client_v4]
name=eduVPN for Linux 4.x (Fedora $releasever)
baseurl=https://app.eduvpn.org/linux/v4/rpm/fedora-$releasever-$basearch
gpgkey=file:///etc/pki/rpm-gpg/RPM-GPG-KEY-python-eduvpn-client_v4
gpgcheck=1
enabled=1
EOF
    sudo dnf install eduvpn-client
    exit 0
}

case $VERSION_CODENAME in
    # ubuntu and deb versions
    "jammy" | "noble" | "plucky" | "questing" | "resolute" | "bookworm" | "trixie")
	install_deb "$VERSION_CODENAME"
	;;
    # redirect linux mint 21.x codenames to jammy
    "vanessa" | "vera" | "victoria")
	install_deb "jammy"
	;;
    # redirect linux mint codenames to noble
    "wilma" | "xia" | "zara" | "zena")
	install_deb "noble"
	;;
    # redirect linux mint 6 codename to bookworm
    "faye")
	install_deb "bookworm"
    ;;
    # redirect linux mint 7 codename to trixie
    "gigi")
	install_deb "trixie"
	;;
esac

# No codename or unsupported codename, get based on name
case $NAME in
    "Fedora Linux")
	install_fedora
	;;
    *)
	echo "OS: \"$NAME\" with codename \"$VERSION_CODENAME\" is not supported"
	exit 1
	;;
esac
