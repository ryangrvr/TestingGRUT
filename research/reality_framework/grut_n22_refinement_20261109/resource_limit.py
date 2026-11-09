"""Bound child committed memory using documented Windows Job Objects."""
import ctypes as C,os,subprocess,time,json
from ctypes import wintypes as W

class Basic(C.Structure):
    _fields_=[('PerProcessUserTimeLimit',C.c_longlong),('PerJobUserTimeLimit',C.c_longlong),('LimitFlags',W.DWORD),
      ('MinimumWorkingSetSize',C.c_size_t),('MaximumWorkingSetSize',C.c_size_t),('ActiveProcessLimit',W.DWORD),
      ('Affinity',C.c_size_t),('PriorityClass',W.DWORD),('SchedulingClass',W.DWORD)]
class IO(C.Structure):
    _fields_=[(n,C.c_ulonglong) for n in ('ReadOperationCount','WriteOperationCount','OtherOperationCount','ReadTransferCount','WriteTransferCount','OtherTransferCount')]
class Extended(C.Structure):
    _fields_=[('BasicLimitInformation',Basic),('IoInfo',IO),('ProcessMemoryLimit',C.c_size_t),('JobMemoryLimit',C.c_size_t),
              ('PeakProcessMemoryUsed',C.c_size_t),('PeakJobMemoryUsed',C.c_size_t)]
class Counters(C.Structure):
    _fields_=[('cb',W.DWORD),('PageFaultCount',W.DWORD)]+[(n,C.c_size_t) for n in
      ('PeakWorkingSetSize','WorkingSetSize','QuotaPeakPagedPoolUsage','QuotaPagedPoolUsage','QuotaPeakNonPagedPoolUsage',
       'QuotaNonPagedPoolUsage','PagefileUsage','PeakPagefileUsage','PrivateUsage')]

class Job:
    def __init__(self,limit):
        if os.name!='nt':raise RuntimeError('This resource-bounded protocol requires Windows')
        self.k=C.WinDLL('kernel32',use_last_error=True);self.p=C.WinDLL('psapi',use_last_error=True)
        self.k.CreateJobObjectW.argtypes=[C.c_void_p,W.LPCWSTR];self.k.CreateJobObjectW.restype=W.HANDLE
        self.k.SetInformationJobObject.argtypes=[W.HANDLE,C.c_int,C.c_void_p,W.DWORD];self.k.SetInformationJobObject.restype=W.BOOL
        self.k.QueryInformationJobObject.argtypes=[W.HANDLE,C.c_int,C.c_void_p,W.DWORD,C.c_void_p];self.k.QueryInformationJobObject.restype=W.BOOL
        self.k.AssignProcessToJobObject.argtypes=[W.HANDLE,W.HANDLE];self.k.AssignProcessToJobObject.restype=W.BOOL
        self.k.CloseHandle.argtypes=[W.HANDLE];self.k.CloseHandle.restype=W.BOOL
        self.p.GetProcessMemoryInfo.argtypes=[W.HANDLE,C.POINTER(Counters),W.DWORD];self.p.GetProcessMemoryInfo.restype=W.BOOL
        self.handle=self.k.CreateJobObjectW(None,None)
        if not self.handle:raise C.WinError(C.get_last_error())
        info=Extended();info.BasicLimitInformation.LimitFlags=0x100|0x2000;info.ProcessMemoryLimit=limit
        if not self.k.SetInformationJobObject(self.handle,9,C.byref(info),C.sizeof(info)):
            self.close();raise C.WinError(C.get_last_error())
    def assign(self,process):
        if not self.k.AssignProcessToJobObject(self.handle,int(process._handle)):raise C.WinError(C.get_last_error())
    def sample(self,process):
        info=Counters();info.cb=C.sizeof(info)
        if not self.p.GetProcessMemoryInfo(int(process._handle),C.byref(info),C.sizeof(info)):return {}
        return {'private_bytes':int(info.PrivateUsage),'working_set_bytes':int(info.WorkingSetSize),'peak_working_set_bytes':int(info.PeakWorkingSetSize)}
    def peak(self):
        info=Extended()
        if not self.k.QueryInformationJobObject(self.handle,9,C.byref(info),C.sizeof(info),None):raise C.WinError(C.get_last_error())
        return int(info.PeakProcessMemoryUsed)
    def close(self):
        if self.handle:self.k.CloseHandle(self.handle);self.handle=None

def bounded_run(cmd,log_path,label,seconds,limits,telemetry_path):
    job=Job(limits['process_memory_bytes']);process=None;start=time.monotonic();last=0.;samples=[];error=None
    try:
        with log_path.open('x',encoding='utf-8') as log:
            process=subprocess.Popen(cmd,stdout=log,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW)
            job.assign(process)
            while process.poll() is None:
                elapsed=time.monotonic()-start
                if elapsed>seconds:error='wall_time_limit';process.kill();break
                sample=job.sample(process)
                if elapsed-last>=limits['heartbeat_seconds']:
                    record={'point':label,'elapsed':elapsed,**sample};samples.append(record)
                    print(json.dumps(record),flush=True);last=elapsed
                time.sleep(limits['poll_seconds'])
            process.wait()
        result={'point':label,'returncode':process.returncode,'wall_seconds':time.monotonic()-start,
          'peak_private_bytes':job.peak(),'error':error,'samples':samples,'memory_limit_bytes':limits['process_memory_bytes']}
    except BaseException as exc:
        if process is not None and process.poll() is None:process.kill();process.wait()
        result={'point':label,'error':str(exc),'wall_seconds':time.monotonic()-start,'samples':samples,'returncode':None}
        raise
    finally:
        job.close()
        if 'result' in locals():
            with telemetry_path.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    return result
