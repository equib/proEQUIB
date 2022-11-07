from ipykernel.kernelbase import Kernel
from IPython.paths import locate_profile
from IPython.core.displaypub import publish_display_data
from pexpect import replwrap,EOF,spawn

import signal
from subprocess import check_output
import tempfile
import re
import os
from glob import glob
from shutil import rmtree
from base64 import b64encode
from distutils.spawn import find_executable

__version__ = '0.5'

version_pat = re.compile(r'Version (\d+(\.\d+)+)')

class GDLKernel(Kernel):
    implementation = 'gdl_kernel'
    implementation_version = __version__
    language = 'GDL'
    @property
    def language_version(self):
        try:
            m = version_pat.search(self.banner)
            return m.group(1)
        except:
            return "Version ?.?"

    _banner = None
    @property
    def banner(self):
        if self._banner is None:
            try:
                self._banner = check_output([self._executable, '--version']).decode('utf-8')
            except:
                self._banner = ''

        return self._banner
    
    language_info = {'name': 'gdl',
                     'codemirror_mode': 'idl',
                     'mimetype': 'text/x-idl',
                     'file_extension': '.pro'}

    def __init__(self, **kwargs):
        Kernel.__init__(self, **kwargs)
        self._start_gdl()

        try:
            self.hist_file = os.path.join(locate_profile(),'gdl_kernel.hist')
        except:
            self.hist_file = None
            self.log.warn('No default profile found, history unavailable')

        self.max_hist_cache = 1000
        self.hist_cache = []

    def _start_gdl(self):
        # Signal handlers are inherited by forked processes, and we can't easily
        # reset it from the subprocess. Since kernelapp ignores SIGINT except in
        # message handlers, we need to temporarily reset the SIGINT handler here
        # so that GDL and its children are interruptible.
        sig = signal.signal(signal.SIGINT, signal.SIG_DFL)
        try:
            self._executable = find_executable("gdl")
            #self._child  = spawn(self._executable+' --use-wx ',timeout = 300, encoding='utf-8')
            self._child  = spawn(self._executable,timeout = 300, encoding='utf-8')
            self.gdlwrapper = replwrap.REPLWrapper(self._child,u"GDL> ",None)
        finally:
            signal.signal(signal.SIGINT, sig)

        self.gdlwrapper.run_command("!quiet=1 & defsysv,'!inline',1 & !more=0".rstrip(), timeout=None)
        # Compile GDL routines/functions
        dirname = os.path.dirname(os.path.abspath(__file__))
        self.gdlwrapper.run_command(".compile "+dirname+"/snapshot.pro",timeout=None)

    def do_execute(self, code, silent, store_history=True, user_expressions=None,
                   allow_stdin=False):

        if not code.strip():
            return {'status': 'ok', 'execution_count': self.execution_count,
                    'payloads': [], 'user_expressions': {}}

        elif (code.strip() == 'exit' or code.strip() == 'quit'):
            self.do_shutdown(False)
            return {'status':'abort','execution_count':self.execution_count}

        elif (code.strip().startswith('.') or code.strip().startswith('@')):
            # This is a GDL Executive command
            output = self.gdlwrapper.run_command(code.strip(), timeout=None) 

            if not silent:
                stream_content = {'name': 'stdout', 'text':output}
                self.send_response(self.iopub_socket, 'stream', stream_content)

            return {'status': 'ok', 'execution_count': self.execution_count,
                    'payloads': [], 'user_expressions': {}}

        if code.strip() and store_history:
            self.hist_cache.append(code.strip())

        interrupted = False
        tfile_code = tempfile.NamedTemporaryFile(mode='w+t',dir=os.path.expanduser("~"))
        tfile_post = tempfile.NamedTemporaryFile(mode='w+t',dir=os.path.expanduser("~"))
        plot_dir = tempfile.mkdtemp(dir=os.path.expanduser("~"))
        plot_format = 'png'

        postcall = """
            if !D.NAME eq 'X' then begin
                device,window_state=winds_arefgij
                if !inline and total(winds_arefgij) ne 0 then begin
                    w_CcjqL6MA = where(winds_arefgij ne 0,nw_CcjqL6MA)
                    for i_KEv8eW6E=0,nw_CcjqL6MA-1 do begin
                        wset, w_CcjqL6MA[i_KEv8eW6E]
                        outfile_c5BXq4dV = '%(plot_dir)s/__fig'+strtrim(i_KEv8eW6E,2)+'.png'
                        ii_rsApk4JS = snapshot(outfile_c5BXq4dV)
                        wdelete
                    endfor
	        endif
	    endif else begin
                if !D.NAME eq 'Z' then begin
                    device
                    if (total(tvrd()) ne 0.0) then begin
                        i_KEv8eW6E=0
                        outfile_c5BXq4dV = '%(plot_dir)s/__fig'+strtrim(i_KEv8eW6E,2)+'.png'
                        ii_rsApk4JS = snapshot(outfile_c5BXq4dV)
                        erase
	            endif
	        endif
	    endelse
        end
        """ % locals()

        try:
            tfile_code.file.write(code.rstrip()+"\na_adfadfw=1\nend")
            tfile_code.file.close()
            tfile_post.file.write(postcall.rstrip())
            tfile_post.file.close()
            output = self.gdlwrapper.run_command(".run "+tfile_code.name, timeout=None)
            self.gdlwrapper.run_command(".run "+tfile_post.name,timeout=None)

            # Publish images if there are any
            images = [open(imgfile, 'rb').read() for imgfile in glob("%s/*.png" % plot_dir)]

            display_data=[]

            for image in images:
                display_data.append({'image/png': b64encode(image).decode('ascii')})

            for data in display_data:
                self.send_response(self.iopub_socket, 'display_data',{'data':data,'metadata':{}})
        except KeyboardInterrupt:
            self.gdlwrapper.child.sendintr()
            interrupted = True
            self.gdlwrapper._expect_prompt()
            output = self.gdlwrapper.child.before
        except EOF:
            output = self.gdlwrapper.child.before + 'Restarting GDL'
            self._start_gdl()
        finally:
            tfile_code.close()
            tfile_post.close()
            rmtree(plot_dir)

        if not silent:
            stream_content = {'name': 'stdout', 'text':output}
            self.send_response(self.iopub_socket, 'stream', stream_content)
        
        if interrupted:
            return {'status': 'abort', 'execution_count': self.execution_count}
        
        try:
            exitcode = int(self.run_command('print,0').rstrip())
        except Exception:
            exitcode = 1

        if exitcode:
            return {'status': 'error', 'execution_count': self.execution_count,
                    'ename': '', 'evalue': str(exitcode), 'traceback': []}
        else:
            return {'status': 'ok', 'execution_count': self.execution_count,
                    'payloads': [], 'user_expressions': {}}

    def do_history(self, hist_access_type, output, raw, session=None,
                   start=None, stop=None, n=None, pattern=None, unique=False):

        if not self.hist_file:
            return {'history': []}

        if not os.path.exists(self.hist_file):
            with open(self.hist_file, 'wb') as f:
                f.write('')

        with open(self.hist_file, 'rb') as f:
            history = f.readlines()

        history = history[:self.max_hist_cache]
        self.hist_cache = history
        self.log.debug('**HISTORY:')
        self.log.debug(history)
        history = [(None, None, h) for h in history]

        return {'history': history}

    def do_shutdown(self, restart):
        self.log.debug("**Shutting down")

        self.gdlwrapper.child.kill(signal.SIGKILL)

        if self.hist_file:
            with open(self.hist_file,'wb') as f:
                data = '\n'.join(self.hist_cache[-self.max_hist_cache:])
                f.write(data.encode('utf-8'))

        return {'status':'ok', 'restart':restart}

if __name__ == '__main__':
    from ipykernel.kernelapp import IPKernelApp
    IPKernelApp.launch_instance(kernel_class=GDLKernel)
