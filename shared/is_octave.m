%% test for GNU octave (vs Matlab)
%from <cite>http://wiki.octave.org/Compatibility</cite>
%
function r = is_octave ()
   persistent x;
   if (isempty (x)) %this test is only done once per matlab/octave invocation
     x = exist ('OCTAVE_VERSION', 'builtin');
   end
   r = x;
 end
