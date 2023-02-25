#! /usr/bin/env python3 -B
import os

def xcode():
    os.makedirs('build/xcode', exist_ok=True)
    os.chdir('build/xcode')
    os.system('cmake ../.. -GXcode -DYOCTO_EMBREE=OFF')
    os.system('open RCM_on_mesh.xcodeproj')

if __name__ == '__main__':
    xcode()
