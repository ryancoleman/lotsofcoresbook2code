.. _making_movies_ase:

=============
Making movies
=============

using recordmydesktop
---------------------

A video tutorial can be produced in the following way:

- change the screen resolution to 1280x1024,

- record the movie of the screen (without sound)
  using recordmydesktop_ (gtk-recordMyDesktop),

- convert the resulting ogv into avi using mencoder::

   mencoder video.ogv -o video.avi -oac copy -ovc lavc

- record and edit the sound track using audacity_:

  - use 44100 Hz for recording and save the final file as *sound.wav*,

  - make sure not to keep the microphone to close to avoid signal peaks,

  - cut microphone signal peaks, insert silences, ...

- edit the movie using avidemux_ (to match the sound track):

  - load the *video.avi*: :menuselection:`File --> Open`, make sure to
    use the following options when editing and saving:
    :menuselection:`Video --> Copy`, :menuselection:`Audio --> Copy`,
    :menuselection:`Format --> AVI`,

  - add the *sound.avi*: :menuselection:`Audio --> Main`
    :menuselection:`Track --> Audio` :menuselection:`Source -->
    External WAV`,

  - cut video frames (or copy and insert still frames to extend the video)
    to match the sound track.
    Set beginning mark and end mark -
    the cut or copy/paste operation applies to the selected region,

  - make sure to save intermediate stages when working on the video
    :menuselection:`File --> Save --> Save Video` (as AVI):

    - avidemux caches the audio track so to match the audio
      to a freshly cut video you can copy the audio file into another name,
      and add the sound track from that name,

    - sometimes when cutting frames avidemux does not allow to set the markers
      correctly, and there is **no** undo the last step in avidemux!

  - save the *video_final.avi* (that matches the sound track),
    and encode it into mpeg4 format using mencoder, using two passes::

     opt="vbitrate=550:mbd=2:dc=10 -vf unsharp=l:0.4:c:0.0:hqdn3d"
     mencoder -ovc lavc -lavcopts vcodec=msmpeg4v2:vpass=1 -nosound -o /dev/null video_final.avi
     mencoder video_final.avi -oac mp3lame -af resample=32000:0:2 -lameopts vbr=3:br=80:mode=3 \
             -ovc lavc -lavcopts acodec=mp3lame:vcodec=msmpeg4v2:vpass=2:$opt \
             -info name="Overview and installation of ASE":artist=CAMd:copyright="CAMd 2009" -o video_final_mpeg4.avi

  - convert *video_final.avi* into a 800x600 *swf* file for streaming::

     ffmpeg -i video_final.avi -pass 1 -s 800x600 -b:a 256k -ar 44100 -ac 1 \
            -vcodec flv -b:v 1200k -g 160 -mbd 2 oi_en_800x600.swf
     ffmpeg -i video_final.avi -pass 2 -s 800x600 -b:a 256k -ar 44100 -ac 1 \
            -vcodec flv -b:v 1200k -g 160 -mbd 2 -y oi_en_800x600.swf

.. _recordmydesktop: http://recordmydesktop.sourceforge.net/
.. _audacity: http://audacity.sourceforge.net/
.. _avidemux: http://www.avidemux.org/

using avconf to collect png files
---------------------------------

Load the trajectory and write the images out as single png files, e. g.:

.. literalinclude:: writepngs.py

In case you do not have avconv, install it (ubuntu)::

  sudo apt-get install libav-tools libavcodec-extra-53 libavdevice-extra-53 libavformat-extra-53 libavutil-extra-51 libpostproc-extra-52 libswscale-extra-2

Convert the png files to a movie (img.mov)::

  avconv -i "%d.png" -r 25 -c:v libx264 -crf 20  -pix_fmt yuv420p img.mov

the options are:

-i "img%d.png" uses these files as the input, %d is a placeholder for the number

-r 25 the desired frame rate, 25 FPS in this case

-c:v libx264 use the h264 codec x264

-crf 20 the video quality, 20 is pretty high, the default is 23

-pix_fmt yuv420p a compatible pixel format
