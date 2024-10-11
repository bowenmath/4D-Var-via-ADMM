% Draw the figures after running the main program

% Show the recovered vorticity field every 10 time steps
for i = 1:NT+1
    if mod(i-1,1) == 0
    omphys = reshape(x(:, i),M,N);
    pcolor(real(omphys)');
    caxis([-omegascale omegascale])
    shading('interp')
    axis equal
    axis off
    title 'Recovered Vorticity'
    colormap(brcol)
    %colorbar('horiz')
    hold on

    hold off
    pause(0.000000001)
    end
end