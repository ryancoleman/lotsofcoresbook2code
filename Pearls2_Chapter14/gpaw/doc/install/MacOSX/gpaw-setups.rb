require 'formula'

class GpawSetups < Formula
  homepage 'https://wiki.fysik.dtu.dk/gpaw/'
  url 'https://wiki.fysik.dtu.dk/gpaw-files/gpaw-setups-0.9.11271.tar.gz'
  sha1 'f9cfbdf50451930cbf7fbd32708f7202e6a7ac7c'

  def install
    Dir.mkdir 'gpaw-setups'
    system 'mv *.gz *.pckl gpaw-setups'
    share.mkpath
    share.install Dir['*']
    ENV.j1  # if your formula's build system can't parallelize
  end

  def test
    system "ls #{share}/gpaw-setups"
  end
end
