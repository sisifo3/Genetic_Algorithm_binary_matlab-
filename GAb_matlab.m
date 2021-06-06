function[X] = rana(~)

close all;
clc;
asd = 1;
%************************************************************************
%set selection_metodo = 1 para usar Roluette Wheel selection y 0 para aleatorio poligamico.*                                                 * 
global selection_metodo;                                               %*
selection_metodo = 0;                                                  %*              
%************************************************************************
%variable global
win_apt = 0;
win_cel = zeros(50,3);
win_cel = num2cell(win_cel);
iter = [];
acum_apt = [];
%En este punto especifico del programa, ya se cuenta con una celda Ce_cm
%con una columna de los numeros binarios la segunda de su valor decimal y 
%una tercera con su valor de aptitud.
[Ce_cm] = poblation(asd);
dec = 0;

while win_apt < .96
    dec = dec + 1;
    [celda_ord] = ord_insertion(Ce_cm);
    [ind_elit] = func_gen_elit(celda_ord);
    [celda_com] = aptitud_d(ind_elit);
    [celda_ord] = ord_insertion(celda_com);
    [aptitud_total] = Glob_aptitud(celda_ord);   
    [win_cel,win_apt] = competition_biology(win_cel,win_apt, aptitud_total,celda_ord);
    Ce_cm = win_cel;
    iter(dec) = dec;
    acum_apt(dec) = win_apt;
    
end
X = win_cel;
xq = 0:pi/32:250*pi;
vq1 = interp1(iter,acum_apt,xq);
plot(iter,acum_apt,'*',xq,vq1,':.');
%contour(iter,acum_apt)















function[Ce_cm] = poblation(~)
arr_crom = randi([0 1], 50, 6);%genera la poblacion en binarios 
xr = bi2de(arr_crom,2);% la pasamos a decimal.

z = zeros(50 ,3);
Ce_cm = num2cell(z); %convertimos el array en celda.
%obtener la aptitud, espera un array de numros decimales. 
[xr1] = aptitud_function(xr);
% es necesario poner los numero binarios en arrays de 6.
n = [];
for lr = 1:50
    for kr = 1:6
        n = [n,arr_crom(lr,kr)];
    end
    
    Ce_cm{lr,1} = n;%agregamos todos los binarios a una array
    n = [];
    Ce_cm{lr,2} = xr(lr);
    Ce_cm{lr,3} = xr1(lr);
end

end


%Sección de funciones necesarias para correr el programa principal
%función de ordenamiento

function[xr1] = aptitud_function(xr)  
%obtenemos el valor de aptitud.
xr1 = zeros(1,50);
for ir = 1:length(xr)
    a = xr(ir);
    if a < 32
    %b = (a/32-1).^(2);   % la funcion fittness para evaluar la aptitud
    b = (1/32 .* (a -32)).^2;
    end
    if a > 31
        b = (1/32 .* (a-31)).^2;
        
    end 
    
    xr1(ir) = b;
   
end
end


function[ind_elit] = func_gen_elit(celda_ord)
za = zeros(50 ,3);
ind_elit = num2cell(za);

for ilr = 1:25
   
    [padre1,padre2] =  sel_pad(celda_ord);
    [hijo1,hijo2] = crz_pa(padre1,padre2);
    [b_1,b_2] = fit_func_elit(padre1,padre2,hijo1,hijo2);
    ind_elit{ilr,1} = b_1;
    ind_elit{(ilr+25),1} = b_2;

end
end

function[celda_ord] = ord_insertion(Ce_cm)

for ls = 1:50
    d = ls;
    while((d > 1) && (Ce_cm{d,3}) > (Ce_cm{d-1,3}))
        %Accomadar aptitud.
        var_temp = Ce_cm{d,3};
        Ce_cm{d,3} = Ce_cm{d-1,3};
        Ce_cm{d-1,3} = var_temp;
       
        %acomadar valor decimial, de acuerdo a aptitud. 
        var_temp1 = Ce_cm{d,2};
        Ce_cm{d,2} = Ce_cm{d-1,2};
        Ce_cm{d-1,2} = var_temp1;
        %acomodar valores binarios, de acuerdo a aptitud. 
        var_temp2 = Ce_cm{d,1};
        Ce_cm{d,1} = Ce_cm{d-1,1};
        Ce_cm{d-1,1} = var_temp2;
        
        d = d-1;
        
    end
end
celda_ord= Ce_cm;
end



function [index] = my_own_RWS(Ce_cmp)
% generamos la probabilidad de que sean seleccionados, esta aumenta 
%dependiendo de su fitness

%creamos un valor prioridad.
vec_prio = zeros(1,50);
for le = 1:50
    vec_prio(le) = Ce_cmp{le,3};
end  


%[1] = previous_probability + (fitness / sum_of_fitness) = 0.0 + (1 / 10) = 0.1
%previous_probability = 0.1
x = vec_prio;
xb = flipud(x); %invertimos los valores 
Xpro = xb/sum(xb);   %Generamos todo la matriz con los resultados (fitness / sum_of_fitness)
%Generamos la probabilidad de ser seleccionado la suma de pre_pro + (fit/sum)
proba = zeros(1,50);
prev_proba = 0;
for km = 1:length(x)
    proba(km) = prev_proba + Xpro(km);
    prev_proba = proba(km);
end
%Escogemos al asar el numero en index que necesitamos.
xbp = flipud(proba);
num_rand = rand;

for ksr = 1:length(xbp)
    if num_rand < xbp(ksr)
        index = ksr;
        return
    end
end

%te da una ubicacion por indice de el numero seleccionado.

   
end

%funcion para seleccionar padres.
function [padre1,padre2] =  sel_pad(Ce_cmp)    
%global selection_metodo 
padre1 = zeros(1,1);
padre2 = zeros(1,1);
for rl = 1:2
    if rl ==1

        if selection_metodo == 1 
            [index] = my_own_RWS(Ce_cmp);
            padre1 = Ce_cmp{index,1};
        else
            [SC] = alt_poligamico(Ce_cmp);
            padre1 = Ce_cmp{SC,1};
        end    
    
    end
    if rl ==2
         
         if selection_metodo == 1 
            [index] = my_own_RWS(Ce_cmp);
            padre2 = Ce_cmp{index,1};
         else
            [SC] = alt_poligamico(Ce_cmp);
            padre2 = Ce_cmp{SC,1};
        end   
                
    end
end
 
end    





%la funcion que genera a los hijos.
% la funcion que hace el cruce de padres
function[hijo1,hijo2] = crz_pa(padre1,padre2)
%es necesario concatenar padre 1 y padre 2;
padres = cat(6,padre1,padre2);

hijo1 =[];
hijo2= [];

for er = 1:3
    hijo1 = [hijo1,padres(er)];
    hijo2 = [hijo2,padres(er+3)];
end

for or = 1:3
    hijo1 = [hijo1,padres(or+9)];
    hijo2 = [hijo2, padres(or+6)];
end    
end        



%en la siguiente funci[on aplicamos elitismo buscando los individuos mas
%aptos.

function [b_1,b_2] = fit_func_elit(padre1,padre2,hijo1,hijo2)

    
%compara aptitud de padres he hijos
%para evaluar el mejor.

%Es necesario convertir en celda, para cuando apliquemos elitismo
%podamos regresar los dos mas aptos en valor binario 
Ce_f = [0,0,0;0,0,0;0,0,0;0,0,0];
Ce_f = num2cell(Ce_f);

%metemos los padres he hijos e forma binaria en la 
%primera columna
Ce_f{1,1}=padre1;
Ce_f{2,1}=padre2;
Ce_f{3,1}=hijo1;
Ce_f{4,1}=hijo2;

%convertir a decimal
padre1 = bi2de(padre1,2);
padre2 = bi2de(padre2,2);
hijo1 = bi2de(hijo1,2); 
hijo2 = bi2de(hijo2,2);
%metemos su valor en decimal en la segunda columna
Ce_f{1,2}=padre1;
Ce_f{2,2}=padre2;
Ce_f{3,2}=hijo1;
Ce_f{4,2}=hijo2;

array = cat(4,padre1,padre2,hijo1,hijo2);

%obtenemos el valor de aptitud.
aptitud_arr = zeros(1,50);
for irs = 1:length(array)
    as = array(irs);
    if as < 32
        bs = (as/32-1).^(2);   % la funcion fittness para evaluar la aptitud
    end
    if as > 31
        bs = (1/32 .* (as-31)).^2;
    end 
    
    
    %as = array(irs);
    %bs = (as/32-1).^(2);   % la funcion fittness para evaluar la aptitud
    aptitud_arr(irs) = bs; 
end

% es necesario aplicar el algoritmo de ordenamiento.
%pero primero es necesario asignar un lugar en la celda

for klr = 1:length(aptitud_arr)
    Ce_f{klr,3} = aptitud_arr(klr);
end

% esta lista la celda con tres columnas, valor en binario, decimal, aptitud
%aplicamos ordenamiento para regresar los dos mas aptos.

for lse = 1:4
    de = lse;
    while((de > 1) && (Ce_f{de,3}) > (Ce_f{de-1,3}))
        %Accomadar aptitud.
        var_tempe = Ce_f{de,3};
        Ce_f{de,3} = Ce_f{de-1,3};
        Ce_f{de-1,3} = var_tempe;
        %d = d-1;
        %acomadar valor decimial, de acuerdo a aptitud. 
        var_temp1e = Ce_f{de,2};
        Ce_f{de,2} = Ce_f{de-1,2};
        Ce_f{de-1,2} = var_temp1e;
        %acomodar valores binarios, de acuerdo a aptitud. 
        var_temp2e = Ce_f{de,1};
        Ce_f{de,1} = Ce_f{de-1,1};
        Ce_f{de-1,1} = var_temp2e;
        
        de = de-1;%disp(aptitud_total)
   
        
    end
end
%nos entrega una variable acomodada por aptitud
%entonces regresamos los dos valores mas aptos.
b_1 = Ce_f{1,1};
b_2 = Ce_f{2,1};
end   
        


%funcion apara convertir de binario a decimal, y dar el valor de aptitud

function[celda_com] = aptitud_d(ind_elit)

    
for tr =1:length(ind_elit)    
    ind_elit{tr,2} = bi2de(ind_elit{tr,1},2);   
end
%obtenemos el valor de aptitud.

for tir = 1:length(ind_elit)
   
    ar = ind_elit{tir,2};
    
    
    %a = xr(ir);
    if ar < 32
        br = (ar/32-1).^(2);   % la funcion fittness para evaluar la aptitud
    end
    
    if ar > 31
        br = (1/32 .* (ar-31)).^2;
        
    end 
    
    ind_elit{tir,3} = br;
    
end
celda_com=ind_elit;

end

%la siguiente funcion es para mutaci[on

function [Celda_mu] = mutation_last (Celda)
%es necesario que este ordenado.
% necesitamos la celda ordenada por aptitud
%despues tomamos el de mejor aptitud y le aplicamos 
%mutacion
for tp = 1:5
    less_apt = zeros(1,5); 
    less_apt(tp) = Celda{(40+tp),1};  
    if less_apt(1)== 0
        less_apt(1) = 1;
    else 
        less_apt(1) = 0;
    end
  
    Celda{(40+tp),1} = less_apt;
    Celda_mu = Celda;

end

end




function[aptitud_total] = Glob_aptitud(celda_com)   
    aptitud_temp = zeros(1,50);
    for fre = 1 : 50
        aptitud_temp(fre) = celda_com{fre,3};
    end    
    sum_aptitud_temp = sum(aptitud_temp);
    aptitud_total = sum_aptitud_temp/50;   
end


function[win_cel,win_apt] = competition_biology(win_cel,win_apt, aptitud_total,celda_ord)
     
     
    if win_apt < aptitud_total
        win_cel = celda_ord;
        win_apt = aptitud_total;
    else
       
        return
    end    
       
end    

function [SC] = alt_poligamico(~)

num_alet  = randi(50);
disp(num_alet)
SC = num_alet;

end

end





