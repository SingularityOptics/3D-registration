%% queue search
% athor ：Bytry
%  time：2016年6月21日07:04:56
% demand：data 
clear all;close all;clc;

% load Left_Hor_unwraped_phase_Standard2_LL.mat
% load Left_Ver_unwraped_phase_Standard2_LL.mat
% load Right_Hor_unwraped_phase_Standard2_LL.mat
% load Right_Ver_unwraped_phase_Standard2_LL.mat
load Left_Hor_unwraped_phase.mat
load Left_Ver_unwraped_phase.mat
load Right_Hor_unwraped_phase.mat
load Right_Ver_unwraped_phase.mat
% load  up1.mat;

unwraped_phase(:,:,1) = Left_Hor_unwraped_phase;
unwraped_phase(:,:,2) = Left_Ver_unwraped_phase;
unwraped_phase(:,:,3) = Right_Hor_unwraped_phase;
unwraped_phase(:,:,4) = Right_Ver_unwraped_phase;
figure(1),imshow(unwraped_phase(:,:,1),[]);    % lhor
figure(2),imshow(unwraped_phase(:,:,2),[]);    %lver
figure(3),imshow(unwraped_phase(:,:,3),[]);    %rhor
figure(4),imshow(unwraped_phase(:,:,4),[]);    %rver
lhor_phase = unwraped_phase(:,:,1);lver_phase = unwraped_phase(:,:,2);
rhor_phase = unwraped_phase(:,:,3);rver_phase = unwraped_phase(:,:,4);

%% find initial row point
disp('点击鼠标左键在图像上选取一个初始测量点：');
figure(1);
[xl,yl,botton] = ginput(1);
hold on; plot(xl(1),yl(1),'b+');
disp('点击鼠标左键在图像上选取测量区域，从上至下、从左至右依次选取4个测量区域边缘点：');
figure(3);
[xr,yr,botton] = ginput(4);
hold on; plot(xr(1),yr(1),'y+');
hold on;plot(xr(2),yr(2),'y+');
hold on;plot(xr(3),yr(3),'y+');
hold on;plot(xr(4),yr(4),'y+');
count1 = 1;
tag1 = 0;
row1 = round(yl(1));
col1 = round( xl(1));
phase_lhor = lhor_phase(row1,col1);
phase_lver = lver_phase(row1,col1);
if (phase_lhor~=0 && phase_lver~=0)
        for row2 = round(yr(1)):round(yr(2))
            for col2 = round(xr(3)):round(xr(4))
                if( rver_phase(row2,col2-1)~=0 && rver_phase(row2,col2)~=0  && rver_phase(row2,col2-1) < phase_lver && rver_phase(row2,col2) > phase_lver)  
                   if(rhor_phase(row2,col2)~=0 && rhor_phase(row2+1,col2)~=0 && rhor_phase(row2,col2)>phase_lhor && rhor_phase(row2+1,col2)<phase_lhor)       
                        ur = col2-1 - (rver_phase(row2,col2-1)-phase_lver)/(rver_phase(row2,col2)-rver_phase(row2,col2-1))*1;            
                        vr = row2+1 + (rhor_phase(row2+1,col2)-phase_lhor)/(rver_phase(row2,col2)-rhor_phase(row2+1,col2));
                        Coord_int(count1,:) = [col1,row1,ur,vr];
                        count1 = count1 + 1;
                        tag1 = 1;
                        break;
                   else if (rhor_phase(row2,col2)~=0 &&  rhor_phase(row2-1,col2)~=0 && rhor_phase(row2-1,col2)>phase_lhor && rhor_phase(row2,col2)<phase_lhor) 
                           ur = col2-1 - (rver_phase(row2,col2-1)-phase_lver)/(rver_phase(row2,col2)-rver_phase(row2,col2-1))*1;
                           vr = row2 + (rhor_phase(row2,col2)-phase_lhor)/(rver_phase(row2-1,col2)-rhor_phase(row2,col2));
                           Coord_int(count1,:) = [col1,row1,ur,vr];
                           count1 = count1 + 1;
                           tag1 = 1;
                           break;
                       end
                   end
                end
            end
            if tag1 == 1
                break;
            end
        end
end
figure(3);hold on; plot(Coord_int(3),Coord_int(4),'r+');

% Coord_int = [582,409,823.087264913529,405.399841715486;]
%% create queue
% create destoryed n-array
center_l_int = [Coord_int(1) Coord_int(2)];
center_r_int = [Coord_int(3) Coord_int(4)];
tag2 = zeros(960,1280); 
queue_head = 1;
queue_tail = 1;
neighbour=[1 0;0 1;-1 0; 0 -1];     %4-connect region
q(:,:,queue_tail) = [round(center_l_int) center_r_int];
queue_tail = queue_tail+1;
[ser1 ser2]=size(neighbour);
count2 = 1;
thr = 30;
tic;
while queue_head~=queue_tail
    center_r_int = q(1, 3:4, queue_head);
    [reg ] = centerserach( center_r_int, thr);
    pix=q(1, 1:2, queue_head);        
    for i = 1:ser1 
        pix1= pix + neighbour(i,:);
        if  (lhor_phase( pix1(2),pix1(1))~=0 && lver_phase(pix1(2),pix1(1))~=0 && tag2(pix1(1),pix1(2)) == 0)  %tag == 0,push queue
            q(1, 1:2,queue_tail) =[pix1(1) pix1(2) ];
            tag2(pix1(1), pix1(2)) = 2;    
            col1 = q(1, 2,queue_tail);
            row1 = q(1, 1,queue_tail);
            phase_lhor = lhor_phase(col1,row1);
            phase_lver = lver_phase(col1,row1);
            [n2 m2 dem] = size(reg);
            for j = 1: dem
                col2 = reg(1,1,j); row2 = reg(1,2,j);
                 if( rver_phase(row2,col2-1)~=0 && rver_phase(row2,col2)~=0 && rver_phase(row2,col2-1) < phase_lver && rver_phase(row2,col2) > phase_lver)        
                           if(rhor_phase(row2,col2)~=0 && rhor_phase(row2+1,col2)~=0 && rhor_phase(row2,col2)>phase_lhor && rhor_phase(row2+1,col2)<phase_lhor)       
                                ur = col2-1 - (rver_phase(row2,col2-1)-phase_lver)/(rver_phase(row2,col2)-rver_phase(row2,col2-1))*1;              %Lagrange interpolation
                                vr = row2+1 + (rhor_phase(row2+1,col2)-phase_lhor)/(rver_phase(row2,col2)-rhor_phase(row2+1,col2));
                                 Coord(count2,:) =  [row1,col1,ur,vr];
                                 count2 = count2 + 1;
                                 q(1, 3:4, queue_tail) = [ur vr];
                                 queue_tail=queue_tail+1;
                               break;
                           else if (rhor_phase(row2,col2)~=0 &&  rhor_phase(row2-1,col2)~=0 && rhor_phase(row2-1,col2)>phase_lhor && rhor_phase(row2,col2)<phase_lhor) 
                                   ur = col2-1 - (rver_phase(row2,col2-1)-phase_lver)/(rver_phase(row2,col2)-rver_phase(row2,col2-1))*1;
                                   vr = row2 + (rhor_phase(row2,col2)-phase_lhor)/(rver_phase(row2-1,col2)-rhor_phase(row2,col2));
                                 Coord(count2,:) = [row1,col1,ur,vr];
                                 count2 = count2 + 1;
                                 q(1, 3:4, queue_tail) = [ur vr];
                                 queue_tail=queue_tail+1;
                                   break;
                               end
                      end
                 end
            end
        end
    end
%remove center
q(:,:,queue_head) = [];
tag2(pix(1),pix(2)) = 1;
queue_tail=queue_tail-1;
end
Coord(count2,:) = Coord_int;
toc;
save test_q.mat Coord;
%% display
figure(1),hold on;
plot(round(Coord(:,1)),round(Coord(:,2)),'r.'); % rhor
figure(3),hold on;
plot(Coord(:,3),Coord(:,4),'g.');  %rver


               
               
                           
                         




