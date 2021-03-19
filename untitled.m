filename   = 'circular_anulus.gif';  % create a new file (attention to not overwrite)    
  del        = 0.1;     
co=figure()
ax1=axes(co)
ax1.XLim=[0 200]
ax1.YLim=[0 200]

for i=20:100;
        M=circular_kernel(i,i-3);
        c=imagesc(ax1,M);xlabel('pixel');ylabel('pixel');title('anular kernel');
       ax1.XLim=[0 200]
       ax1.YLim=[0 200]
       drawnow 
       frame      = getframe(ax1);
       im         = frame2im(frame);
       [imind,cm] = rgb2ind(im,256);   % This creates a matrix of index and color map (cm)
       if i == 20;
         imwrite(imind,cm,filename,'gif','Loopcount',inf,'DelayTime',del);
       else
         imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',del);
       end      
 end