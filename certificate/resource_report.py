"""Portable resource reporting; not used by any mathematical predicate."""
import resource
import sys

def peak_rss_report(platform_name=None, raw=None):
    name = sys.platform if platform_name is None else platform_name
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss if raw is None else raw
    unit = 'bytes' if name == 'darwin' else 'KiB' if name.startswith('linux') else 'platform-native-unknown'
    byte_value = int(value) if name == 'darwin' else int(value) * 1024 if name.startswith('linux') else None
    return {'peak_rss_bytes': byte_value, 'ru_maxrss_raw': value,
            'ru_maxrss_native_unit': unit, 'resource_platform': name}
