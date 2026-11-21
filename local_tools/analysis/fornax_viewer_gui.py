"""Desktop launcher for the local/offline Fornax browser viewer."""
from configparser import Error as ConfigError, RawConfigParser
from datetime import date
import json
from pathlib import Path
from queue import Empty, Queue
import re
import subprocess
import sys
import threading
import webbrowser

from visualize_simulation import LABELS, SNAPSHOT


SCRIPT = Path(__file__).resolve().with_name('visualize_simulation.py')
MODEL = re.compile(r'(Fornax|MW|Sculptor)(\d+)(.*)')


def workbench():
    for root in (*SCRIPT.parents, Path.cwd(), *Path.cwd().parents):
        path = root/'workflow.json'
        if path.is_file() and not path.is_symlink() and (root/'analysis').is_dir():
            return root
    return Path.cwd()


def model_dirs(root):
    found = [(MODEL.fullmatch(p.name), p) for p in root.iterdir() if p.is_dir()]
    return [p for m,p in sorted(((m,p) for m,p in found if m),
                                key=lambda a:(a[0][1], int(a[0][2]), a[0][3]))]


def candidate_inputs(model):
    """Scan only local model data/standard branches, not release or cache trees."""
    directories = [model, model/'isolation', model/'output', *sorted(model.glob('orbit_*'))]
    for directory in tuple(directories):
        if directory.is_dir():
            directories.extend(sorted(directory.glob('snapdir_*')))
            if directory.name.startswith('orbit_'):
                directories.append(directory/'output')
    found = {}
    for directory in directories:
        if not directory.is_dir():
            continue
        for p in sorted(directory.iterdir()):
            if not p.is_file():
                continue
            match = SNAPSHOT.fullmatch(p.name)
            if match:
                if match[3] not in (None, '0'):
                    continue
                key = (1, int(match[2]), str(p))
            elif p.suffix.lower()=='.ini':
                with p.open(errors='replace') as stream:
                    head = stream.readline(256).split()
                if len(head)!=4 or not all(s.isdigit() for s in head):
                    continue
                key = (0, 0, str(p))
            else:
                continue
            found.setdefault(p.resolve(),(key, p))
    return [path for key,path in sorted(found.values())]


def source_info(path):
    lines = [str(path)]
    if path.suffix.lower()=='.ini':
        with path.open() as stream:
            count = stream.readline(256).split()
        if len(count)!=4 or not all(s.isdigit() for s in count):
            raise ValueError('不是粒子初态 INI（不能使用参数配方）')
        lines.append('ASCII 粒子初态；气体 / TP1 / 旧星 / TP3：'+ ' / '.join(count))
        for p in path.parent.glob('*build.json'):
            record = json.loads(p.read_text())
            named = record.get('local_file_prefix','') + str(record.get('particle_ic',''))
            if named and (path.parent/named).resolve() == path.resolve():
                lines.append(f'构建：{record.get("build_id") or "历史记录"}；seed={record.get("seed")}')
                break
    else:
        import h5py
        with h5py.File(path,'r') as f:
            a = f['Header'].attrs
            lines.append(f'GIZMO HDF5；代码时间 {a.get("Time","—")}；分片数 {a.get("NumFilesPerSnapshot",1)}；'
                         f'HubbleParam={a.get("HubbleParam","—")}（不自动换算）')
            counts = a.get('NumPart_Total',a.get('NumPart_ThisFile'))
            high = a.get('NumPart_Total_HighWord', [0]*6)
            if counts is not None:
                lines.append('；'.join(f'{label}: {int(n)+(int(h)<<32):,}' for label,n,h in zip(LABELS,counts,high)
                                      if int(n) or int(h)) or '头信息未记录非零粒子计数')
    return '\n'.join(lines)


def model_info(model):
    """Give a recipe summary even if the model has no saved local particles."""
    recipe = model.parent/'IniFile'/f'IC_{model.name}.ini'
    if not recipe.is_file():
        return f'{model.name}：没有登记的本地配方。'
    p = RawConfigParser()
    p.read(recipe)
    items = []
    for suffix, label in (('stellar','恒星'),('gas','气体')):
        d, c = f'Density_{suffix}', f'Component_{suffix}'
        if p.has_section(d):
            values=[]
            mass=p.get(d,'mass',fallback=None)
            if mass is not None:
                values.append(f'配置质量 {float(mass)*1e10:.4g} M☉')
            count=p.get(c,'nbody',fallback=None)
            if count is not None:
                values.append(f'{int(count):,} 粒子')
            scale=p.get(d,'scaleRadius',fallback=None)
            if scale is not None:
                values.append(f'尺度 {scale} kpc')
            if values:
                items.append(label+'：'+ ' / '.join(values))
    return model.name+' 配方（不是快照实测） · '+'；'.join(items) if items else str(recipe)


def launch():
    import tkinter as tk
    from tkinter import filedialog, messagebox, ttk

    root = workbench()
    try:
        app = tk.Tk()
    except tk.TclError as error:
        raise RuntimeError(f'无法打开桌面窗口：{error}') from error
    app.title('Fornax 模拟查看器')
    app.geometry('970x730')
    app.minsize(830,650)
    content = ttk.Frame(app, padding=14)
    content.pack(fill='both', expand=True)
    content.columnconfigure(0, weight=1)
    content.rowconfigure(2, weight=1)
    bar = ttk.Frame(content)
    bar.grid(row=0,column=0,sticky='ew')
    bar.columnconfigure(1,weight=1)
    ttk.Label(bar,text='模型').grid(row=0,column=0,padx=(0,9))
    models = model_dirs(root)
    current = tk.StringVar()
    choose = ttk.Combobox(bar,textvariable=current,values=[p.name for p in models])
    choose.grid(row=0,column=1,sticky='ew')
    summary = tk.StringVar(value='Ctrl/Shift 可多选；默认只列出本机可读的快照和粒子初态。')
    ttk.Label(content,textvariable=summary,wraplength=910).grid(row=1,column=0,sticky='w',pady=8)
    files = ttk.Frame(content)
    files.grid(row=2,column=0,sticky='nsew')
    files.columnconfigure(0,weight=1)
    files.rowconfigure(0,weight=1)
    listing = tk.Listbox(files,selectmode='extended',exportselection=False,height=11)
    listing.grid(row=0,column=0,sticky='nsew')
    scroll = ttk.Scrollbar(files,command=listing.yview)
    scroll.grid(row=0,column=1,sticky='ns')
    listing.configure(yscrollcommand=scroll.set)
    info = tk.StringVar(value='请选择模型，或直接选择本地文件。')
    ttk.Label(content,textvariable=info,wraplength=910,justify='left').grid(row=3,column=0,sticky='w',pady=9)
    opts = ttk.Frame(content)
    opts.grid(row=4,column=0,sticky='ew')
    for i in (1,3,5):
        opts.columnconfigure(i,weight=1)
    radius, axis, units = [tk.StringVar(value=v) for v in ('5','short','fornax')]
    points, light, ml, types = [tk.StringVar(value=v) for v in ('2500','20000','2.6','0,2,4')]
    center, bulk = tk.StringVar(),tk.StringVar()

    def opt(label,var,row,col,choices=None):
        ttk.Label(opts,text=label).grid(row=row,column=col,sticky='w',padx=(0,5),pady=3)
        field = ttk.Combobox(opts,textvariable=var,values=choices,state='readonly',width=18) if choices else ttk.Entry(opts,textvariable=var,width=17)
        field.grid(row=row,column=col+1,sticky='ew',padx=(0,13),pady=3)

    opt('球形孔径',radius,0,0); opt('参考轴',axis,0,2,('short','native-z','angular-momentum')); opt('单位',units,0,4,('fornax','native'))
    opt('每组显示点数',points,1,0); opt('恒星光采样数',light,1,2); opt('V 波段 M/L',ml,1,4)
    opt('固定中心 XYZ',center,2,0); opt('固定速度 XYZ',bulk,2,2); opt('显示粒子类型',types,2,4)
    ttk.Label(content,text='留空自动用 TP2 定中心；正式轨道需核对中心。面亮度预览不等同于观测拟合。').grid(row=5,column=0,sticky='w',pady=6)
    actions = ttk.Frame(content)
    actions.grid(row=6,column=0,sticky='ew')
    status = tk.StringVar(value='就绪')
    ttk.Label(actions,textvariable=status).pack(side='left')
    candidates = []
    done = Queue()
    last = None

    def update_info(_=None):
        try:
            sel = listing.curselection()
            if not sel:
                return
            desc = source_info(candidates[sel[-1]])
            if len(sel)>1:
                desc += f'\n已选 {len(sel)} 帧；下面的头信息对应最后选中的文件。'
            info.set(desc)
        except (OSError,ValueError,KeyError,ImportError) as error:
            info.set(str(error))

    def populate(paths):
        candidates[:] = paths
        listing.delete(0,'end')
        for p in paths:
            listing.insert('end',str(p.relative_to(root)) if p.is_relative_to(root) else str(p))
        if paths:
            listing.selection_set(len(paths)-1)
            listing.see('end')
            update_info()
        else:
            info.set('此模型没有找到可读的本地 HDF5 快照或粒子初态 INI。可手动选文件。')

    def update_models(_=None):
        text = current.get().strip()
        choose['values']=[m.name for m in models if text.lower() in m.name.lower()]

    def pick_model(_=None):
        text=current.get().strip()
        if text.isdigit():
            text='Fornax'+text
        selected=[p for p in models if p.name.lower()==text.lower()]
        if selected:
            try:
                current.set(selected[0].name)
                try:
                    summary.set(model_info(selected[0]))
                except (ConfigError,OSError,ValueError) as error:
                    summary.set(f'无法读取配方摘要：{error}')
                populate(candidate_inputs(selected[0]))
            except (ValueError,OSError) as error:
                info.set(str(error))
        else:
            info.set('请选择完整模型名，例如 Fornax8004。')

    def pick_file():
        files = filedialog.askopenfilenames(title='选本地快照或初态',
                     filetypes=[('GIZMO 快照或粒子 INI','*.hdf5 *.h5 *.ini'),('所有文件','*')])
        if files:
            summary.set('手动选取的本地数据；以下头信息来自实际文件。')
            populate([Path(f) for f in files])
            listing.selection_set(0,'end')
            update_info()

    def run_job(command,out):
        try:
            run = subprocess.run(command,cwd=root,text=True,capture_output=True,check=False)
            done.put((run.returncode==0,run.stdout if run.returncode==0 else run.stderr,out))
        except OSError as error:
            done.put((False,str(error),out))

    def finish():
        nonlocal last
        try:
            success,message,out = done.get_nowait()
        except Empty:
            app.after(150,finish)
            return
        view.configure(state='normal')
        save.configure(state='normal')
        if success:
            last = out
            status.set(message.strip())
            again.configure(state='normal')
            webbrowser.open(out.as_uri())
        else:
            status.set('生成失败')
            messagebox.showerror('无法生成',message.strip())

    def output_path(chosen):
        name=re.sub(r'[^A-Za-z0-9_.-]','_',current.get() or chosen[0].stem)
        folder = root/'analysis'/f'{date.today().isoformat()}_3D_{name}'
        default = folder/'viewer.html'
        if default.exists():
            i=2
            while (folder/f'viewer-{i}.html').exists():
                i+=1
            default=folder/f'viewer-{i}.html'
        return default

    def start(custom=False):
        chosen = [candidates[i] for i in listing.curselection()]
        if not chosen:
            status.set('请至少选择一个文件。')
            return
        out=output_path(chosen)
        if custom:
            selected=filedialog.asksaveasfilename(title='保存离线交互网页',initialdir=str(out.parent if out.parent.is_dir() else root/'analysis'),
                    initialfile=out.name,defaultextension='.html',filetypes=[('HTML','*.html')])
            if not selected:
                return
            out=Path(selected).resolve()
        try:
            aperture=float(radius.get())
        except ValueError:
            status.set('孔径必须是数字。')
            return
        command=[sys.executable,str(SCRIPT),*map(str,chosen),'--output',str(out),'--radius',radius.get(),
                 '--axis-radius',str(min(2.,aperture)),'--axis',axis.get(),'--units',units.get(),'--points',points.get(),
                 '--profile-points',light.get(),'--mass-to-light',ml.get(),'--types',types.get()]
        for key,value in (('--center',center.get()),('--bulk',bulk.get())):
            if value.strip():
                command.extend([key,*re.split(r'[\s,]+',value.strip())])
        view.configure(state='disabled')
        save.configure(state='disabled')
        status.set('正在读取数据并生成网页…')
        threading.Thread(target=run_job,args=(command,out),daemon=True).start()
        app.after(150,finish)

    ttk.Button(bar,text='手动选文件',command=pick_file).grid(row=0,column=2,padx=9)
    view=ttk.Button(actions,text='生成并打开',command=start)
    view.pack(side='right')
    save=ttk.Button(actions,text='另存为…',command=lambda:start(True))
    save.pack(side='right',padx=8)
    again=ttk.Button(actions,text='打开上一次结果',state='disabled',command=lambda:webbrowser.open(last.as_uri()))
    again.pack(side='right',padx=8)
    listing.bind('<<ListboxSelect>>',update_info)
    choose.bind('<KeyRelease>',update_models)
    choose.bind('<<ComboboxSelected>>',pick_model)
    choose.bind('<Return>',pick_model)
    choose.focus_set()
    app.mainloop()
    return 0
