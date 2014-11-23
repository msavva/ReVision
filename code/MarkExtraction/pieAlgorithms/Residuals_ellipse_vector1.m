function [RSS, XYproj] = Residuals_ellipse_vector1(XY,ParG)
%
%   Projecting a given set of points onto an ellipse
%   and computing the distances from the points to the ellipse
%
%   This is a modified version of an iterative algorithm published by D. Eberly
%     Internet publication: "Distance from a point to an ellipse in 2D" (2004)
%                           Geometric Tools, LLC, www.geometrictools.com
%     Book publication: "3D Game Engine Design", 2nd edition.
%                       Morgan Kaufmann Publishers, San Francisco, CA, 2007.
%                              (see Section 14.13.1)
%
%   Input:  XY(n,2) is the array of coordinates of n points x(i)=XY(i,1), y(i)=XY(i,2)
%           ParG is a vector 5x1 of the ellipse parameters
%           ParG =  [Center(1:2), Axes(1:2), Angle]'
%                Center - the coordinates of the ellipse's center
%                Axes   - the axes (major, minor)
%                Angle  - the angle of tilt of the ellipse
%
%   Output:  RSS is the Residual Sum of Squares (the sum of squares of the distances)
%            XYproj is the array of coordinates of projections
%
%   The algorithm is proven to converge and reaches an accuracy of 7-8 significant digit
%   It takes 4-5 iterations per point, on average.


Center = ParG(1:2);   Axes = ParG(3:4);  Angle = ParG(5);

n = size(XY,1);
XYproj = zeros(n,2);
tolerance = 1e-9;

%  First handling the circle case

if abs((Axes(1)-Axes(2))/Axes(1))<tolerance
    phiall = angle(XY(:,1)-Center(1) + sqrt(-1)*(XY(:,2)-Center(2)));
    XYproj = [Axes(1)*cos(phiall)+Center(1) Axes(2)*sin(phiall)+Center(2)];
    RSS = sqrt(sum((XY-XYproj).^2,2));
    return;
end

%  Now dealing with proper ellipses

a = Axes(1); b = Axes(2);
aa = a^2;  bb = b^2;

tol_a = tolerance*a;
tol_b = tolerance*b;
tol_aa = tolerance*aa;

%  Matrix Q for rotating the points and the ellipse to the canonical system

s = sin(Angle); c = cos(Angle);
Q = [c -s; s c];

%  data points in canonical coordinates

XY0  = [XY(:,1)-Center(1) XY(:,2)-Center(2)]*Q;
XYA = abs(XY0);
Tini = max(a*(XYA(:,1)-a),b*(XYA(:,2)-b));

%  main loop over the data points

z1 = zeros(n,1);
z2 = zeros(n,2);

z1(XYA(:,1) == 0) = 1;
z1(XYA(:,1) ~= 0) = sign(XY0(XYA(:,1) ~= 0));
z2(XYA(:,2) == 0) = 1;
z2(XYA(:,2) ~= 0) = sign(XY0(XYA(:,2) ~= 0));

u = XYA(:,1); v = XYA(:,2);
ua = u*a; vb = v*b;

u_lt_tola = u < tol_a;
XYproj(u_lt_tola & XY0(:,2) < 0,:) = repmat([0,-b], [length(find(u_lt_tola & XY0(:,2) < 0)) 1]);
XYproj(u_lt_tola & XY0(:,2) >= 0,:) = repmat([0,b], [length(find(u_lt_tola & XY0(:,2) >= 0)) 1]);

v_lt_tolb = v < tol_b;
u_lt_abb = u < a-bb/a;
xproj = (aa/(aa-bb))*u;
xproj_a = xproj(v_lt_tolb & u_lt_abb);
z1_a = z1(v_lt_tolb & u_lt_abb);
z2_a = z2(v_lt_tolb & u_lt_abb);
z1_b = z1(v_lt_tolb & ~u_lt_abb);
% XYproj(v_lt_tolb & u_lt_abb, :) = repmat([z1_a.*xproj_a z2_a.*b.*sqrt(1-(xproj_a./a).^2)], ...
%                                   [length(find(v_lt_tolb & u_lt_abb)),1]);
% XYproj(v_lt_tolb & ~u_lt_abb, :) = repmat([z1_b.*a repmat(0, [length(find(v_lt_tolb & ~u_lt_abb)),1])], ...
%                                   [length(find(v_lt_tolb & ~u_lt_abb)),1]);
XYproj(v_lt_tolb & u_lt_abb, :) = [z1_a.*xproj_a z2_a.*b.*sqrt(1-(xproj_a./a).^2)];
XYproj(v_lt_tolb & ~u_lt_abb, :) = [z1_b.*a repmat(0, [length(find(v_lt_tolb & ~u_lt_abb)),1])];                            

other_points = ~u_lt_tola & ~v_lt_tolb;
T = Tini(other_points);
for iter=1:100
    Taa = T + aa;
    Tbb = T + bb;
    PP1 = (ua(other_points)./Taa).^2;
    PP2 = (vb(other_points)./Tbb).^2;
    F  = PP1 + PP2 - 1;
    if F<0; break; end;
    Fder = 2*(PP1./Taa + PP2./Tbb);
    Ratio = F./Fder;
    if (Ratio<tol_aa); break; end;
    T = T + Ratio;
end
xproj = (XY0(other_points,1)*aa)./Taa;
XYproj(other_points,:) = [xproj sign(XY0(other_points,2)).*b.*sqrt(1-(xproj./a).^2)];

% for i=1:n
%     u = XYA(i,1);  v = XYA(i,2);
%     ua = u*a;      vb = v*b;
% 
% %     if (u == 0)
% %         z1 = 1;
% %     else
% %         z1 = sign(XY0(i,1));
% %     end
% % 
% %     if (v == 0)
% %         z2 = 1;
% %     else
% %         z2 = sign(XY0(i,2));
% %     end
% 
%     %       does the point lie on the minor axis?
% 
% %     if u<tol_a
% %         if XY0(i,2)<0
% %             XYproj(i,:) = [0 -b];
% %         else
% %             XYproj(i,:) = [0 b];
% %         end
% %         continue;
% %     end
% 
%     %       does the point lie on the major axis?
% 
% %     if v<tol_b
% %         if u<a-bb/a
% %             xproj = aa*u/(aa-bb);
% %             XYproj(i,:) = [z1*xproj z2*b*sqrt(1-(xproj/a)^2)];
% %         else
% %             XYproj(i,:) = [z1*a 0];
% %         end
% %         continue;
% %     end
% 
%     %      generic case: start the iterative procedure
% 
%     T = Tini(i);
%     for iter=1:100
%         Taa = T + aa;
%         Tbb = T + bb;
%         PP1 = (ua/Taa)^2;
%         PP2 = (vb/Tbb)^2;
%         F  = PP1 + PP2 - 1;
%         if F<0; break; end;
%         Fder = 2*(PP1/Taa + PP2/Tbb);
%         Ratio = F/Fder;
%         if (Ratio<tol_aa); break; end;
%         T = T + Ratio;
%     end
% 
%     %       compute the projection of the point onto the ellipse
% 
%     xproj = XY0(i,1)*aa/Taa;
%     XYproj(i,:) = [xproj sign(XY0(i,2))*b*sqrt(1-(xproj/a)^2)];
% 
% end % end the main loop

%    rotate back to the original system
XYproj = XYproj*Q';
XYproj = [XYproj(:,1)+Center(1) XYproj(:,2)+Center(2)];
RSS = sqrt(sum((XY-XYproj).^2,2));
%RSS = norm(XY-XYproj,'fro')^2;
end   % Residuals_ellipse
