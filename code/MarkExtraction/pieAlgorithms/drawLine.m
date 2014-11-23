function drawLine(line_params, color)
    if ~exist('color','var')
        color = 'green';
    end
    hold on;
    xy = [line_params.x0; line_params.x0+line_params.d];
    plot(xy(:,1),xy(:,2),'LineWidth',2,'Color',color);
    plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
    plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');
    hold off;
end