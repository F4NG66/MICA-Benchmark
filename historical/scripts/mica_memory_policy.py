"""Respect both host RAM availability and the effective container memory limit."""
from pathlib import Path

def memory_snapshot(meminfo=Path('/proc/meminfo'),cgroup=Path('/sys/fs/cgroup')):
 values=dict(line.split(':',1) for line in meminfo.read_text().splitlines())
 host=int(values['MemAvailable'].split()[0])*1024
 snapshot={'host_available_bytes':host,'effective_available_bytes':host}
 limit=cgroup/'memory.max';current=cgroup/'memory.current'
 if limit.exists() and current.exists():
  value=limit.read_text().strip()
  if value!='max':
   cap=int(value);used=int(current.read_text().strip());remaining=max(0,cap-used)
   snapshot.update(cgroup_limit_bytes=cap,cgroup_current_bytes=used,cgroup_headroom_bytes=remaining,effective_available_bytes=min(host,remaining))
 return snapshot

def available_memory():return memory_snapshot()['effective_available_bytes']
