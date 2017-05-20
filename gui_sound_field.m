function varargout = gui_sound_field(varargin)
% GUI_SOUND_FIELD MATLAB code for gui_sound_field.fig
%      GUI_SOUND_FIELD, by itself, creates a new GUI_SOUND_FIELD or raises the existing
%      singleton*.
%
%      H = GUI_SOUND_FIELD returns the handle to a new GUI_SOUND_FIELD or the handle to
%      the existing singleton*.
%
%      GUI_SOUND_FIELD('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_SOUND_FIELD.M with the given input arguments.
%
%      GUI_SOUND_FIELD('Property','Value',...) creates a new GUI_SOUND_FIELD or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_sound_field_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_sound_field_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_sound_field

% Last Modified by GUIDE v2.5 27-Oct-2015 14:38:02

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_sound_field_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_sound_field_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before gui_sound_field is made visible.
function gui_sound_field_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_sound_field (see VARARGIN)
% Choose default command line output for gui_sound_field
set(handles.quadrate_checkbox,'value',0);
set(handles.circular_checkbox,'value',1);

%background_image1 = importdata('background.jpg');
%axes(handles.background_axes);
%image(background_image1);
%alpha(0.5)
%axis off
sign1 = imread('tubiao.png');
axes(handles.sign_axes);
image(sign1);
axis off
plot_image=importdata('plot.png');
set(handles.plot_pushbutton,'CDATA',plot_image)
play_image=importdata('stop2.png');
set(handles.play_pushbutton,'CDATA',play_image)
set(handles.text23,'String',['dB/m']);
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gui_sound_field wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = gui_sound_field_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in plot_pushbutton.
function plot_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to plot_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

quadrate_flag=get(handles.quadrate_checkbox,'value');
circular_flag=get(handles.circular_checkbox,'value');
z0_double=str2num(get(handles.z0_edit,'string'));
f0=str2num(get(handles.frequency_edit,'string'));                  %  Transducer center frequency [MHz] 
f0=f0*1e6;
c=str2num(get(handles.c_edit,'string'));                           %  Speed of sound [m/s]( 水 1500 或者 钢 5900) 
decay=str2num(get(handles.decay_edit,'string'));                   %衰减系数
%------判断输入参数的正误--------%
if f0<20000
    errorcall= errordlg( 'The Ultrasonic Frequency is Higher than 20000Hz at least. Please input the parameter again.' , 'Error'  ) ;
    dbcont;
end
if c<300
    errorcall= errordlg( 'The Speed of Sound is too small. Please input the parameter again.' , 'Error'  ) ;
    dbcont;
end
if decay<0
    errorcall= errordlg( 'The Attenuation of Sound may not be less than zero. Please input the parameter again.' , 'Error'  ) ;
    dbcont;
end
%------判断输入参数的正误--------%
%-----clear axes-----%
axes(handles.acoustic_axis_axes);
acoustic_axis_delete=get(gca,'children');
delete(acoustic_axis_delete);
axes(handles.sound_field_axes);
sound_field_delete=get(gca,'children');
delete(sound_field_delete);
axes(handles.directivity_axes);
directivity_delete=get(gca,'children');
delete(directivity_delete);
%-------------------%


pausetime=0.05;
flag=1;
p0=101.325*10^3;                             %大气压强101.325*10^3pa
Rp=1;
lambda=c/f0;                      %  Wavelength 
k=2*pi/lambda;                    %波数
w=2*pi*f0;                        %角频率
if quadrate_flag==0 && circular_flag==1
    %---------------------------------------圆形活塞----------------------------------------------%
    %圆形活塞参数
    Rs=str2num(get(handles.radius_edit,'string'));                      %活塞半径25mm
    Rs=Rs*10^-3;
    Fs=pi*Rs^2;                       %活塞的面积
    %各种材料的声阻抗
    Z_air=340*1.29;
    Z_water=1500*10^3;
    Z_steel=str2num(get(handles.Z_R_edit,'string')); 
    %工件参数
    z_d=str2num(get(handles.z_d_edit,'string'));
    z_d=z_d*10^-3;                     %厚度
    %定义网格，数值求解
    n=150;
    x=linspace(-1.5*Rs,1.5*Rs,n);
    y=linspace(-1.5*Rs,1.5*Rs,n);
    %轴向初始距离
    near_field=(2*Rs)^2/(4*lambda);                                        %近场距离
    z0=z0_double*near_field;
    %------判断输入参数的正误--------%
    if Rs<=0
    errorcall= errordlg( 'The Size of Instrument may not be less than zero. Please input the parameter again.' , 'Error'  ) ;
    dbcont;
    end
    if z_d<=0
    errorcall= errordlg( 'Workpiece Size may not be less than zero. Please input the parameter again.' , 'Error'  ) ;
    dbcont;
    end
    if z0_double<0 || z0>z_d
    errorcall= errordlg( 'The Initial Position may not be less than zero or more than Workpiece Size. Please input the parameter again.' , 'Error'  ) ;
    dbcont;
    end
    %------判断输入参数的正误--------%
    %反馈近场距离到界面
    set(handles.near_field_text,'String',num2str(near_field*10^3));
    %--------------------------%
    z=z0;
    detaz=z;
    [x,y]=meshgrid(x,y);
    r=sqrt(x.^2+y.^2+z.^2);
    Sita=acos(z./r);
    Phi=asin(y./(r.*sin(Sita)));
    R_acoustic_axis_incident=linspace(0,z_d,10000);
    R_acoustic_axis_reflect=linspace(z_d,10*(2*Rs)^2/(4*lambda),10000);
    %-----轴线上声压-----%
    P_field_incident=P_circular_acoustic_axis(lambda,Rs,R_acoustic_axis_incident).*exp(-decay*R_acoustic_axis_incident/8.68);
    P_field_reflect=rp(Z_steel,Z_air)*P_circular_acoustic_axis(lambda,Rs,R_acoustic_axis_reflect).*exp(-decay*R_acoustic_axis_reflect/8.68);
    %由结论，直接代入参数
    P=Rp*Sound_pressure_circular( k,Rs,Sita, w,p0,c,r).*exp(-decay*r/8.68);
    P=abs(P);
    %指向性
    [theta,phi]=meshgrid(linspace(0,2*pi,3*n),linspace(0,pi/2,3*n));
    X=k*Rs*sin(phi);
    J1=besselj(1,X);
    D=abs(2*J1./X);
    %坐标变换
    %[Dx,Dy,Dz]=[D.*((sin(Phi)).*cos(Sita)),D.*((sin(Sita)).*sin(Phi)),D.*((cos(Phi)).*cos(Sita-Sita))];
    Dx=D.*((sin(phi)).*cos(theta));
    Dy=D.*((sin(theta)).*sin(phi));
    Dz=D.*((cos(phi)).*cos(theta-theta));

    axes(handles.acoustic_axis_axes);
    plot(R_acoustic_axis_incident,abs(P_field_incident));
    hold on
    plot(R_acoustic_axis_reflect,abs(P_field_reflect),'r');
    plot(linspace(z_d,z_d,10),linspace(0,2,10),'--ks','LineWidth',1,'MarkerSize',2);
    text(z_d,abs(P_field_incident(1000)),['reflect point']);
    grid on;
    title(['Circular pressure transducer acoustic axis of distribution'],'FontSize',7);
    xlabel(['z/m']);
    ylabel(['P/P0']);

    axes(handles.sound_field_axes);
    p=surf(x,y,P+0.7*((1-flag)*z_d+flag*z)*10^8);
    set(p,'facealpha',z0/z-0.1);
    material shiny
    shading interp
    colormap(cool)
    light ('position',[-1 -0.5 2],'style','infinite')
    hold on
    cube=surf_cube( 4*Rs,4*Rs,0.7*z_d*10^8,-2*Rs,-2*Rs,0);
    axis([-3*Rs 3*Rs -3*Rs 3*Rs 0 1.5*(0.7*z_d*10^8)]);
    title(['When z=',num2str(((1-flag)*z_d+flag*z)*10^3),'mm,Circular pressure transducer acoustic axis of the three-dimensional distribution'],'FontSize',7);
    text(0,0,0.04*z_d*10^8,['Workpiece:Plane of incidence'],'color','r');
    text(0,0,0.7*z_d*10^8,['Reflective surface'],'color','r');

    axes(handles.directivity_axes);
    D_surf=surf(Dx,Dy,Dz);
    set(D_surf,'facealpha',0.5);
    axis([-0.2 0.2 -0.2 0.2 0 1])
    title(['Beam directiviy'],'FontSize',7);
    colorbar
    material shiny
    shading interp
    colormap(cool)
    %---------------------------------------圆形活塞----------------------------------------------%
elseif quadrate_flag==1 && circular_flag==0
    %---------------------------------------矩形活塞----------------------------------------------%   
    %各种材料的声阻抗
    Z_air=340*1.29;
    Z_water=1500*10^3;
    Z_steel=str2num(get(handles.Z_R_edit,'string'));

    %矩形活塞换能器参数:边长a b，面积Fs
    a=str2num(get(handles.length_edit,'string'))*10^-3;
    b=str2num(get(handles.width_edit,'string'))*10^-3;
    Fs=a*b;
    %工件参数
    z_d=str2num(get(handles.z_d_edit,'string'));
    z_d=z_d*10^-3;                     %厚度
     %---------近场-------------%
    near_field=Fs/(2*pi*lambda);
    z0=z0_double*Fs/near_field;
    %------判断输入参数的正误--------%
    if a<=0 || b<=0
    errorcall= errordlg( 'The Size of Instrument may not be less than zero. Please input the parameter again.' , 'Error'  ) ;
    dbcont;
    end
    if z_d<=0
    errorcall= errordlg( 'Workpiece Size may not be less than zero. Please input the parameter again.' , 'Error'  ) ;
    dbcont;
    end
    if z0_double<0 || z0>z_d
    errorcall= errordlg( 'The Initial Position may not be less than zero or more than Workpiece Size. Please input the parameter again.' , 'Error'  ) ;
    dbcont;
    end
    %------判断输入参数的正误--------%
    %空间网格化求数值解
    n=120;                            %
    x=linspace(-a,a,n);
    y=linspace(-b,b,n);
    Sita=linspace(0,pi/2,n);
    apha=linspace(0,2*pi,n);
    [Sita,apha]=meshgrid(Sita,apha);
    
    X_acoustic_axis=linspace(-a,a,n);
    Z_acoustic_axis_incident=linspace(0,z_d,n);
    Z_acoustic_axis_reflect=linspace(z_d,1.5*z_d,n);
    [X_acoustic_axis_incident,Z_acoustic_axis_incident]=meshgrid(X_acoustic_axis,Z_acoustic_axis_incident);
    [X_acoustic_axis_reflect,Z_acoustic_axis_reflect]=meshgrid(X_acoustic_axis,Z_acoustic_axis_reflect);
    R_acoustic_axis_incident=sqrt(X_acoustic_axis_incident.^2+Z_acoustic_axis_incident.^2);
    R_acoustic_axis_reflect=sqrt(X_acoustic_axis_reflect.^2+Z_acoustic_axis_reflect.^2);
    phi_acoustic_axis=0;
    sita_acoustic_axis_incident=asin(X_acoustic_axis_incident./(R_acoustic_axis_incident.*cos(phi_acoustic_axis)));
    sita_acoustic_axis_reflect=asin(X_acoustic_axis_reflect./(R_acoustic_axis_reflect.*cos(phi_acoustic_axis)));
    %-------------------%
    %-----------反馈近场距离到界面--------------%
    set(handles.near_field_text,'String',num2str(near_field*10^3));
    %------------------------------------------%
    z=z0;
    [x,y]=meshgrid(x,y);
    r=sqrt(x.^2+y.^2+z.^2);
    phi=asin(y./r);
    sita=asin(x./(r.*cos(phi)));
    A=k*a*cos(apha).*sin(Sita);
    B=k*b*sin(apha).*sin(Sita);
    %直接用结论求解
    %轴线声压
    P_acoustic_axis_incident=P_rectangle_acoustic_axis( Fs,a,k,lambda,R_acoustic_axis_incident ,sita_acoustic_axis_incident,phi_acoustic_axis);
    P_acoustic_axis_incident=P_acoustic_axis_incident.*exp(-decay*R_acoustic_axis_incident/8.68);
    P_acoustic_axis_reflect=P_rectangle_acoustic_axis( Fs,a,k,lambda,R_acoustic_axis_reflect ,sita_acoustic_axis_reflect,phi_acoustic_axis);
    P_acoustic_axis_reflect=rp(Z_steel,Z_air)*P_acoustic_axis_reflect.*exp(-decay*R_acoustic_axis_reflect/8.68);
    axes(handles.acoustic_axis_axes);
    p_acoustic_axis_inciden=surf(X_acoustic_axis_incident,Z_acoustic_axis_incident,abs(P_acoustic_axis_incident));
    hold on
    p_acoustic_axis_reflect=surf(X_acoustic_axis_reflect,Z_acoustic_axis_reflect,abs(P_acoustic_axis_reflect));
    plot3(linspace(0,0,10),linspace(z_d,z_d,10),linspace(0,20,10),'--ks','LineWidth',1,'MarkerSize',2);
    text(0,z_d,15,['reflect point']);
    title(['Recrangular pressure transducer acoustic axis of distribution'],'FontSize',7);
    material shiny
    shading interp
    colormap(cool)
    set(p_acoustic_axis_inciden,'facealpha',0.7);
    set(p_acoustic_axis_reflect,'facealpha',0.7);
    light ('position',[1 1 20],'style','infinite');
    view([80,16]);
    %声压
    P=Sound_pressure_rectangle( p0,Fs,lambda,r,k,a,b,sita,phi).*exp(-decay*r/8.68);
    P=abs(P);
    %指向性
    D_aph_Sita_omiga=(sin(A)./A).*(sin(B)./B);
    D_aph_Sita_omiga=abs(D_aph_Sita_omiga);
    axes(handles.sound_field_axes);
    p=surf(x,y,P+0.7*((1-flag)*z_d+flag*z)*10^7);
    %set(p,'facealpha',z0/z-0.2);
    title(['When z=',num2str(((1-flag)*z_d+flag*z)*10^3),'mm,Recrangular pressure transducer acoustic axis of the three-dimensional distribution'],'FontSize',7);
    material shiny
    shading interp
    colormap(cool)
    light ('position',[-1 -0.5 2],'style','infinite')
    hold on
    cube=surf_cube( 4*a,4*b,0.7*z_d*10^7,-2*a,-2*b,0);
    axis([-3*a 3*a -3*b 3*b 0 1.1*0.7*z_d*10^7]);
    %坐标变换
    Dx=D_aph_Sita_omiga.*(sin(Sita)).*cos(apha);
    Dy=D_aph_Sita_omiga.*(sin(Sita)).*sin(apha);
    Dz=D_aph_Sita_omiga.*(cos(Sita)).*cos(apha-apha);
    text(0,0,0,['Workpiece:Plane of incidence'],'color','r');
    text(0,0,0.7*z_d*10^7,['Reflective surface'],'color','r');
    
    axes(handles.directivity_axes);
    D_surf=surf(Dx,Dy,Dz);
    set(D_surf,'facealpha',0.5)
    axis([-0.2 0.2 -0.2 0.2 0 1])
    title(['Beam directiviy'],'FontSize',7);
    material shiny
    shading interp
    colormap(cool)
    colorbar
end


% --- Executes on button press in play_pushbutton.
function play_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to play_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
quadrate_flag=get(handles.quadrate_checkbox,'value');
circular_flag=get(handles.circular_checkbox,'value');
z0_double=str2num(get(handles.z0_edit,'string'));
f0=str2num(get(handles.frequency_edit,'string'));                  %  Transducer center frequency [MHz] 
f0=f0*1e6;
c=str2num(get(handles.c_edit,'string'));                           %  Speed of sound [m/s]( 水 1500 或者 钢 5900) 
decay=str2num(get(handles.decay_edit,'string'));                   %衰减系数

%------判断输入参数的正误--------%
if f0<20000
    errorcall= errordlg( 'The Ultrasonic Frequency is Higher than 20000Hz at least. Please input the parameter again.' , 'Error'  ) ;
    dbcont;
end
if c<300
    errorcall= errordlg( 'The Speed of Sound is too small. Please input the parameter again.' , 'Error'  ) ;
    dbcont;
end
if decay<0
    errorcall= errordlg( 'The Attenuation of Sound may not be less than zero. Please input the parameter again.' , 'Error'  ) ;
    dbcont;
end
%------判断输入参数的正误--------%
%-----clear axes-----%
axes(handles.acoustic_axis_axes);
acoustic_axis_delete=get(gca,'children');
delete(acoustic_axis_delete);
axes(handles.sound_field_axes);
sound_field_delete=get(gca,'children');
delete(sound_field_delete);
axes(handles.directivity_axes);
directivity_delete=get(gca,'children');
delete(directivity_delete);
%-------------------%


pausetime=0.05;
flag=1;

p0=101.325*10^3;                             %大气压强101.325*10^3pa
Rp=1;
lambda=c/f0;                      %  Wavelength 
k=2*pi/lambda;                    %波数
w=2*pi*f0;                        %角频率

if quadrate_flag==0 && circular_flag==1
    %---------------------------------------圆形活塞----------------------------------------------%
    %圆形活塞参数
    Rs=str2num(get(handles.radius_edit,'string'));                      %活塞半径25mm
    Rs=Rs*10^-3;
    Fs=pi*Rs^2;                       %活塞的面积
    %各种材料的声阻抗
    Z_air=340*1.29;
    Z_water=1500*10^3;
    Z_steel=str2num(get(handles.Z_R_edit,'string'));
    %工件参数
    z_d=str2num(get(handles.z_d_edit,'string'));
    z_d=z_d*10^-3;                     %厚度
    %定义网格，数值求解
    n=150;
    x=linspace(-1.5*Rs,1.5*Rs,n);
    y=linspace(-1.5*Rs,1.5*Rs,n);
    %轴向初始距离
    near_field=(2*Rs)^2/(4*lambda);                                        %近场距离
    z0=z0_double*near_field;
    %------判断输入参数的正误--------%
    if Rs<=0
    errorcall= errordlg( 'The Size of Instrument may not be less than zero. Please input the parameter again.' , 'Error'  ) ;
    dbcont;
    end
    if z_d<=0
    errorcall= errordlg( 'Workpiece Size may not be less than zero. Please input the parameter again.' , 'Error'  ) ;
    dbcont;
    end
    if z0_double<0 || z0>z_d
    errorcall= errordlg( 'The Initial Position may not be less than zero or more than Workpiece Size. Please input the parameter again.' , 'Error'  ) ;
    dbcont;
    end
    %------判断输入参数的正误--------%
    %反馈近场距离到界面
    set(handles.near_field_text,'String',num2str(near_field*10^3));
    %--------------------------
    z=z0;
    detaz=z;
    [x,y]=meshgrid(x,y);
    r=sqrt(x.^2+y.^2+z.^2);
    Sita=acos(z./r);
    Phi=asin(y./(r.*sin(Sita)));
    R_acoustic_axis_incident=linspace(0,z_d,10000);
    R_acoustic_axis_reflect=linspace(z_d,10*(2*Rs)^2/(4*lambda),10000);
    %-----轴线上声压-----%
    P_field_incident=P_circular_acoustic_axis(lambda,Rs,R_acoustic_axis_incident).*exp(-decay*R_acoustic_axis_incident/8.68);
    P_field_reflect=rp(Z_steel,Z_air)*P_circular_acoustic_axis(lambda,Rs,R_acoustic_axis_reflect).*exp(-decay*R_acoustic_axis_reflect/8.68);
    %由结论，直接代入参数
    P=Sound_pressure_circular( k,Rs,Sita, w,p0,c,r).*exp(-decay*r/8.68);
    P=abs(P);
    %指向性
    [theta,phi]=meshgrid(linspace(0,2*pi,3*n),linspace(0,pi/2,3*n));
    X=k*Rs*sin(phi);
    J1=besselj(1,X);
    D=abs(2*J1./X);
    %坐标变换
    %[Dx,Dy,Dz]=[D.*((sin(Phi)).*cos(Sita)),D.*((sin(Sita)).*sin(Phi)),D.*((cos(Phi)).*cos(Sita-Sita))];
    Dx=D.*((sin(phi)).*cos(theta));
    Dy=D.*((sin(theta)).*sin(phi));
    Dz=D.*((cos(phi)).*cos(theta-theta));

    axes(handles.acoustic_axis_axes);
    plot(R_acoustic_axis_incident,abs(P_field_incident));
    hold on
    plot(R_acoustic_axis_reflect,abs(P_field_reflect),'r');
    plot(linspace(z_d,z_d,10),linspace(0,2,10),'--ks','LineWidth',1,'MarkerSize',2);
    text(z_d,abs(P_field_incident(1000)),['reflect point']);
    grid on;
    title(['Circular pressure transducer acoustic axis of distribution'],'FontSize',7);
    xlabel(['z/m']);
    ylabel(['P/P0']);

    axes(handles.sound_field_axes);
    p=surf(x,y,P+0.7*((1-flag)*z_d+flag*z)*10^8);
    set(p,'facealpha',z0/z-0.1);
    material shiny
    shading interp
    colormap(cool)
    light ('position',[-1 -0.5 2],'style','infinite')
    hold on
    cube=surf_cube( 4*Rs,4*Rs,0.7*z_d*10^8,-2*Rs,-2*Rs,0);
    axis([-3*Rs 3*Rs -3*Rs 3*Rs 0 1.5*(0.7*z_d*10^8)]);
    text(0,0,0.04*z_d*10^8,['Workpiece:Plane of incidence'],'color','r');
    text(0,0,0.7*z_d*10^8,['Reflective surface'],'color','r');

    axes(handles.directivity_axes);
    D_surf=surf(Dx,Dy,Dz);
    set(D_surf,'facealpha',0.5);
    axis([-0.2 0.2 -0.2 0.2 0 1])
    title(['Beam directiviy'],'FontSize',7);
    colorbar
    material shiny
    shading interp
    colormap(cool)


    while 1
        axes(handles.sound_field_axes);
        r=sqrt(x.^2+y.^2+z.^2);
        Sita=acos(z./r);
        Phi=asin(y./(r.*sin(Sita)));
    %----判断入射波和回波声场模拟---%
    if z>=z_d && flag==1
        %回波声场
        %反射
        Rp=rp(Z_steel,Z_air);
        %作为回波的初始声压值
        P0=P;
        flag=-1;
    end
    z=z+0.005*pausetime;
    %入(反)射波声场
    P=Rp*Sound_pressure_circular( k,Rs,Sita, w,p0,c,r).*exp(-decay*r/8.68);
    P=flag*abs(P);
        %axis([-1.5*a 1.5*a -1.5*b 1.5*b 0 5*p0]);
    %--------------------%
    set(p,'xdata',x,'ydata',y,'zdata',P+0.7*((1-flag)*z_d+flag*z)*10^8)
    set(p,'facealpha',(z0/z)^0.1);
    axis([-3*Rs 3*Rs -3*Rs 3*Rs 0 1.5*(0.7*z_d*10^8)]);

    if z>=2*z_d
        break;
    end
    title(['When z=',num2str(((1-flag)*z_d+flag*z)*10^3),'mm,Circular pressure transducer acoustic axis of the three-dimensional distribution'],'FontSize',7);
    xlabel('x/m');
    ylabel('y/m');
    zlabel('p/pa or z/mm');
    pause(pausetime);
    drawnow
    end
    title(['When z=',num2str(((1-flag)*z_d+flag*z)*10^3),'mm,Circular pressure transducer acoustic axis of the three-dimensional distribution'],'FontSize',7);
    xlabel('x/m');
    ylabel('y/m');
    zlabel('p/pa or z/mm');
    %---------------------------------------圆形活塞----------------------------------------------%
elseif quadrate_flag==1 && circular_flag==0
    %---------------------------------------矩形活塞----------------------------------------------%   
    %各种材料的声阻抗
    Z_air=340*1.29;
    Z_water=1500*10^3;
    Z_steel=str2num(get(handles.Z_R_edit,'string'));

    %矩形活塞换能器参数:边长a b，面积Fs
    a=str2num(get(handles.length_edit,'string'))*10^-3;
    b=str2num(get(handles.width_edit,'string'))*10^-3;
    Fs=a*b;
    %工件参数
    z_d=str2num(get(handles.z_d_edit,'string'));
    z_d=z_d*10^-3;                     %厚度
    %---------近场-------------%
    near_field=Fs/(2*pi*lambda);
    z0=z0_double*Fs/near_field;
    %------判断输入参数的正误--------%
    if a<=0 || b<=0
    errorcall= errordlg( 'The Size of Instrument may not be less than zero. Please input the parameter again.' , 'Error'  ) ;
    dbcont;
    end
    if z_d<=0
    errorcall= errordlg( 'Workpiece Size may not be less than zero. Please input the parameter again.' , 'Error'  ) ;
    dbcont;
    end
    if z0_double<0 || z0>z_d
    errorcall= errordlg( 'The Initial Position may not be less than zero or more than Workpiece Size. Please input the parameter again.' , 'Error'  ) ;
    dbcont;
    end
    %------判断输入参数的正误--------%
    %空间网格化求数值解
    n=120;                            %
    x=linspace(-a,a,n);
    y=linspace(-b,b,n);
    Sita=linspace(0,pi/2,n);
    apha=linspace(0,2*pi,n);
    [Sita,apha]=meshgrid(Sita,apha);
    
   X_acoustic_axis=linspace(-a,a,n);
    Z_acoustic_axis_incident=linspace(0,z_d,n);
    Z_acoustic_axis_reflect=linspace(z_d,1.5*z_d,n);
    [X_acoustic_axis_incident,Z_acoustic_axis_incident]=meshgrid(X_acoustic_axis,Z_acoustic_axis_incident);
    [X_acoustic_axis_reflect,Z_acoustic_axis_reflect]=meshgrid(X_acoustic_axis,Z_acoustic_axis_reflect);
    R_acoustic_axis_incident=sqrt(X_acoustic_axis_incident.^2+Z_acoustic_axis_incident.^2);
    R_acoustic_axis_reflect=sqrt(X_acoustic_axis_reflect.^2+Z_acoustic_axis_reflect.^2);
    phi_acoustic_axis=0;
    sita_acoustic_axis_incident=asin(X_acoustic_axis_incident./(R_acoustic_axis_incident.*cos(phi_acoustic_axis)));
    sita_acoustic_axis_reflect=asin(X_acoustic_axis_reflect./(R_acoustic_axis_reflect.*cos(phi_acoustic_axis)));
    %-------------------%
    %-----------反馈近场距离到界面--------------%
    set(handles.near_field_text,'String',num2str(near_field*10^3));
    %------------------------------------------%
    z=z0;
    [x,y]=meshgrid(x,y);
    r=sqrt(x.^2+y.^2+z.^2);
    phi=asin(y./r);
    sita=asin(x./(r.*cos(phi)));
    A=k*a*cos(apha).*sin(Sita);
    B=k*b*sin(apha).*sin(Sita);
    %直接用结论求解
    %轴线声压
    P_acoustic_axis_incident=P_rectangle_acoustic_axis( Fs,a,k,lambda,R_acoustic_axis_incident ,sita_acoustic_axis_incident,phi_acoustic_axis);
    P_acoustic_axis_incident=P_acoustic_axis_incident.*exp(-decay*R_acoustic_axis_incident/8.68);
    P_acoustic_axis_reflect=P_rectangle_acoustic_axis( Fs,a,k,lambda,R_acoustic_axis_reflect ,sita_acoustic_axis_reflect,phi_acoustic_axis);
    P_acoustic_axis_reflect=rp(Z_steel,Z_air)*P_acoustic_axis_reflect.*exp(-decay*R_acoustic_axis_reflect/8.68);
    axes(handles.acoustic_axis_axes);
    p_acoustic_axis_inciden=surf(X_acoustic_axis_incident,Z_acoustic_axis_incident,abs(P_acoustic_axis_incident));
    hold on
    p_acoustic_axis_reflect=surf(X_acoustic_axis_reflect,Z_acoustic_axis_reflect,abs(P_acoustic_axis_reflect));
    plot3(linspace(0,0,10),linspace(z_d,z_d,10),linspace(0,20,10),'--ks','LineWidth',1,'MarkerSize',2);
    text(0,z_d,15,['reflect point']);
    title(['Recrangular pressure transducer acoustic axis of distribution'],'FontSize',7);
    material shiny
    shading interp
    colormap(cool)
    set(p_acoustic_axis_inciden,'facealpha',0.7);
    set(p_acoustic_axis_reflect,'facealpha',0.7);
    light ('position',[1 1 20],'style','infinite');
    view([80,16]);
    %声压
    P=Sound_pressure_rectangle( p0,Fs,lambda,r,k,a,b,sita,phi).*exp(-decay*r/8.68);
    P=abs(P);
    %指向性
    D_aph_Sita_omiga=(sin(A)./A).*(sin(B)./B);
    D_aph_Sita_omiga=abs(D_aph_Sita_omiga);
    axes(handles.sound_field_axes);
    p=surf(x,y,P+0.7*((1-flag)*z_d+flag*z)*10^7);
    %set(p,'facealpha',z0/z-0.2);
    material shiny
    shading interp
    colormap(cool)
    light ('position',[-1 -0.5 2],'style','infinite')
    hold on
    cube=surf_cube( 4*a,4*b,0.7*z_d*10^7,-2*a,-2*b,0);
    %坐标变换
    Dx=D_aph_Sita_omiga.*(sin(Sita)).*cos(apha);
    Dy=D_aph_Sita_omiga.*(sin(Sita)).*sin(apha);
    Dz=D_aph_Sita_omiga.*(cos(Sita)).*cos(apha-apha);
    text(0,0,0,['Workpiece:Plane of incidence'],'color','r');
    text(0,0,0.7*z_d*10^7,['Reflective surface'],'color','r');
    
    
    axes(handles.directivity_axes);
    D_surf=surf(Dx,Dy,Dz);
    set(D_surf,'facealpha',0.5)
    axis([-0.2 0.2 -0.2 0.2 0 1])
    title(['Beam directiviy'],'FontSize',7);
    material shiny
    shading interp
    colormap(cool)
    colorbar
    
    while 1
    axes(handles.sound_field_axes);
    r=sqrt(x.^2+y.^2+z.^2);
    phi=asin(y./r);
    sita=asin(x./(r.*cos(phi)));
    %----判断入射波和回波声场模拟---%
    if z>=z_d && flag==1
        %回波声场
        %反射
        Rp=rp(Z_steel,Z_air);
        flag=-1;
    end
    z=z+0.005*pausetime;
    P=Rp*Sound_pressure_rectangle( p0,Fs,lambda,r,k,a,b,sita,phi).*exp(-decay*r/8.68);
    P=flag*abs(P);
        %axis([-1.5*a 1.5*a -1.5*b 1.5*b 0 5*p0]);

    %--------------------%
    set(p,'xdata',x,'ydata',y,'zdata',P+0.7*((1-flag)*z_d+flag*z)*10^7)
    set(p,'facealpha',(z0/z)^0.1);
    axis([-3*a 3*a -3*b 3*b 0 1.1*0.7*z_d*10^7]);

    if z>=2*z_d
        break;
    end
    title(['When z=',num2str(((1-flag)*z_d+flag*z)*10^3),'mm,Recrangular pressure transducer acoustic axis of the three-dimensional distribution'],'FontSize',7)
    xlabel('x/m');
    ylabel('y/m');
    zlabel('p/pa or z/mm');
    pause(pausetime);
    drawnow
    end
    title(['When z=',num2str(((1-flag)*z_d+flag*z)*10^3),'mm,Recrangular pressure transducer acoustic axis of the three-dimensional distribution'],'FontSize',7)
    xlabel('x/m');
    ylabel('y/m');
    zlabel('p/pa or z/mm');

    %---------------------------------------矩形活塞----------------------------------------------%   

    
end




function z_d_edit_Callback(hObject, eventdata, handles)
% hObject    handle to z_d_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of z_d_edit as text
%        str2double(get(hObject,'String')) returns contents of z_d_edit as a double


% --- Executes during object creation, after setting all properties.
function z_d_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to z_d_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Z_R_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Z_R_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Z_R_edit as text
%        str2double(get(hObject,'String')) returns contents of Z_R_edit as a double


% --- Executes during object creation, after setting all properties.
function Z_R_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Z_R_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function decay_edit_Callback(hObject, eventdata, handles)
% hObject    handle to decay_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of decay_edit as text
%        str2double(get(hObject,'String')) returns contents of decay_edit as a double


% --- Executes during object creation, after setting all properties.
function decay_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to decay_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function c_edit_Callback(hObject, eventdata, handles)
% hObject    handle to c_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of c_edit as text
%        str2double(get(hObject,'String')) returns contents of c_edit as a double


% --- Executes during object creation, after setting all properties.
function c_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to c_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in quadrate_checkbox.
function quadrate_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to quadrate_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.circular_checkbox,'value',0);

% Hint: get(hObject,'Value') returns toggle state of quadrate_checkbox


% --- Executes on button press in circular_checkbox.
function circular_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to circular_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.quadrate_checkbox,'value',0);

% Hint: get(hObject,'Value') returns toggle state of circular_checkbox



function frequency_edit_Callback(hObject, eventdata, handles)
% hObject    handle to frequency_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of frequency_edit as text
%        str2double(get(hObject,'String')) returns contents of frequency_edit as a double


% --- Executes during object creation, after setting all properties.
function frequency_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frequency_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function radius_edit_Callback(hObject, eventdata, handles)
% hObject    handle to radius_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of radius_edit as text
%        str2double(get(hObject,'String')) returns contents of radius_edit as a double


% --- Executes during object creation, after setting all properties.
function radius_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to radius_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function length_edit_Callback(hObject, eventdata, handles)
% hObject    handle to length_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of length_edit as text
%        str2double(get(hObject,'String')) returns contents of length_edit as a double


% --- Executes during object creation, after setting all properties.
function length_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to length_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function width_edit_Callback(hObject, eventdata, handles)
% hObject    handle to width_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of width_edit as text
%        str2double(get(hObject,'String')) returns contents of width_edit as a double


% --- Executes during object creation, after setting all properties.
function width_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to width_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function file_menu_Callback(hObject, eventdata, handles)
% hObject    handle to file_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function sound_field_menu_Callback(hObject, eventdata, handles)
% hObject    handle to sound_field_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function close_menu_Callback(hObject, eventdata, handles)
% hObject    handle to close_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function z0_edit_Callback(hObject, eventdata, handles)
% hObject    handle to z0_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of z0_edit as text
%        str2double(get(hObject,'String')) returns contents of z0_edit as a double


% --- Executes during object creation, after setting all properties.
function z0_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to z0_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
