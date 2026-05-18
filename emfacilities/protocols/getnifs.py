#!/usr/bin/python
# Python module for listing network interfaces with name,
# index and addresses
# Based on getifaddrs.py from pydlnadms [http://code.google.com/p/pydlnadms/].
# Only tested on Linux!

import sys
from socket import AF_INET, AF_INET6, inet_ntop
from ctypes import (
    Structure, Union, POINTER,
    pointer, get_errno, cast,
    c_ushort, c_byte, c_void_p, c_char_p, c_uint, c_uint16, c_uint32
)
import ctypes.util
import ctypes


class struct_sockaddr(Structure):
    _fields_ = [
        ('sa_family', c_ushort),
        ('sa_data', c_byte * 14), ]


class struct_sockaddr_in(Structure):
    _fields_ = [
        ('sin_family', c_ushort),
        ('sin_port', c_uint16),
        ('sin_addr', c_byte * 4)]


class struct_sockaddr_in6(Structure):
    _fields_ = [
        ('sin6_family', c_ushort),
        ('sin6_port', c_uint16),
        ('sin6_flowinfo', c_uint32),
        ('sin6_addr', c_byte * 16),
        ('sin6_scope_id', c_uint32)]


class union_ifa_ifu(Union):
    _fields_ = [
        ('ifu_broadaddr', POINTER(struct_sockaddr)),
        ('ifu_dstaddr', POINTER(struct_sockaddr)), ]


class struct_ifaddrs(Structure):
    pass


struct_ifaddrs._fields_ = [
    ('ifa_next', POINTER(struct_ifaddrs)),
    ('ifa_name', c_char_p),
    ('ifa_flags', c_uint),
    ('ifa_addr', POINTER(struct_sockaddr)),
    ('ifa_netmask', POINTER(struct_sockaddr)),
    ('ifa_ifu', union_ifa_ifu),
    ('ifa_data', c_void_p), ]

libc = ctypes.CDLL(ctypes.util.find_library('c'))


def ifap_iter(ifap):
    ifa = ifap.contents
    while True:
        yield ifa
        if not ifa.ifa_next:
            break
        ifa = ifa.ifa_next.contents


def getfamaddr(sa):
    family = sa.sa_family
    addr = None
    if family == AF_INET:
        sa = cast(pointer(sa), POINTER(struct_sockaddr_in)).contents
        addr = inet_ntop(family, sa.sin_addr)
    elif family == AF_INET6:
        sa = cast(pointer(sa), POINTER(struct_sockaddr_in6)).contents
        addr = inet_ntop(family, sa.sin6_addr)
    return family, addr


class NetworkInterface(object):
    """
    Represents low level network interface information used to identify and describe
    the communication endpoints available on a system. The collection of related
    structures and helper utilities provides access to IPv4 and IPv6 interface
    addresses, interface identifiers, and address family information through native
    operating system networking services.

    AI Generated:

    Network Interface Discovery (struct_sockaddr) - User Manual
        Overview

        This utility provides a platform-oriented mechanism for discovering and
        inspecting the network interfaces available on a machine. Its primary goal
        is to expose interface names, numerical identifiers, and associated network
        addresses in a format that can be consumed by higher level applications,
        monitoring services, or distributed computing environments.

        In practical environments, network interface inspection is important when
        configuring communication between services, validating connectivity, or
        selecting the correct interface for data transfer. Systems with multiple
        adapters, virtual interfaces, VPN connections, or containerized networking
        often require explicit interface discovery to ensure that communication
        occurs through the intended route.

        General Workflow

        The utility interacts with the operating system networking layer to retrieve
        the complete list of active interfaces currently registered on the machine.
        Each detected interface is represented with a human readable name together
        with its associated addressing information. Both IPv4 and IPv6 environments
        are supported, allowing compatibility with modern dual-stack infrastructures.

        The interface discovery process is designed to provide a lightweight snapshot
        of the networking state at execution time. This makes it suitable for
        initialization tasks, runtime diagnostics, and adaptive network selection in
        distributed applications.

        IPv4 and IPv6 Support

        Modern computing environments commonly expose both IPv4 and IPv6 addresses.
        The utility separates and interprets these address families independently so
        that applications can determine which protocol versions are available on each
        interface.

        IPv4 addresses are generally used in traditional local area networks and
        legacy infrastructures, while IPv6 addresses are increasingly important in
        cloud environments, large institutional networks, and modern internet-facing
        deployments. Supporting both standards ensures broader compatibility and
        simplifies deployment across heterogeneous systems.

        Interface Identification

        Each interface is associated with a numerical index in addition to its name.
        These identifiers are especially useful in low level networking workflows,
        multicast communication, routing configuration, and socket-based applications
        where interfaces must be referenced unambiguously.

        In systems with many virtual or dynamically created interfaces, relying only
        on textual names may not always be sufficient. Numerical identifiers provide
        a stable mechanism for programmatic interaction with the networking stack.

        Error Handling and Robustness

        The utility is designed to tolerate partially inaccessible or invalid
        interfaces while continuing discovery of the remaining network configuration.
        This behavior is valuable in environments where interfaces may appear or
        disappear dynamically, such as cloud orchestration systems, container
        platforms, or hardware undergoing reconfiguration.

        When an interface cannot be interpreted correctly, the discovery process
        continues without interrupting the overall inspection workflow. This improves
        robustness in heterogeneous operating system environments.

        Outputs and Interpretation

        The resulting interface objects provide a concise representation of the
        network configuration available on the host. Each object contains the
        interface name, its system index, and the available IP addresses grouped by
        protocol family.

        These outputs can be used directly for diagnostics, connection management,
        distributed service discovery, or automated infrastructure validation. In
        many workflows, they also serve as the basis for selecting the preferred
        communication interface in multi-network systems.

        Practical Recommendations

        In production systems, it is often useful to verify that the expected IPv4
        or IPv6 addresses are present before starting network-dependent services.
        Interfaces associated with virtual machines, VPNs, or container overlays
        should be inspected carefully because they may introduce additional routing
        paths that influence communication behavior.

        When working in distributed computing or scientific processing environments,
        selecting the correct interface can significantly improve reliability and
        prevent accidental exposure of services through unintended networks.

        Final Perspective

        Network interface discovery is a foundational capability for many networking
        and distributed computing tasks. By exposing structured access to interface
        names, protocol families, and address information, this utility simplifies
        communication setup and helps applications adapt to complex and changing
        networking environments.
    """


def __init__(self, name):
    self.name = name
    self.index = libc.if_nametoindex(name)
    self.addresses = {}


def __str__(self):
    return "%s [index=%d, IPv4=%s, IPv6=%s]" % (
        self.name, self.index,
        self.addresses.get(AF_INET),
        self.addresses.get(AF_INET6))


def getName(self):
    return self.name


def getIndex(self):
    return self.index


def getAddresses(self):
    return self.addresses


def get_network_interfaces():
    ifap = POINTER(struct_ifaddrs)()
    result = libc.getifaddrs(pointer(ifap))
    if result != 0:
        raise OSError(get_errno())
    del result
    retval = {}
    for ifa in ifap_iter(ifap):
        try:
            name = ifa.ifa_name.decode("utf-8")
            i = retval.get(name)
            if not i:
                i = retval[name] = NetworkInterface(name)
                family, addr = getfamaddr(ifa.ifa_addr.contents)
                if addr:
                    i.addresses[family] = addr
        except ValueError:
            del retval[name]
            print("get_network_interfaces: Can not connect to NIC %s" % name)
        except Exception as ex:
            print("Unexpected error:", ex)
    # print retval.values()
    libc.freeifaddrs(ifap)
    return retval.values()


if __name__ == '__main__':
    print([str(ni) for ni in get_network_interfaces()])
