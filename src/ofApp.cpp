#include <iostream>
#include "ofApp.h"

//--------------------------------------------------------------
void ofApp::setup()
{
    std::cout << "begin setup" << std::endl;

    Graph initial_graph = Graph::example_graph_AxB(11, 7);
    circle_packing = CirclePacking(initial_graph);

    circle_packing.set_fonts("SpaceMono-Regular.ttf", "SpaceMono-Italic.ttf", "SpaceMono-Bold.ttf", 16);

    ofSetCircleResolution(100);

    std::cout << "end   setup" << std::endl;
}

//--------------------------------------------------------------
void ofApp::update() {}

//--------------------------------------------------------------
void ofApp::draw()
{
    ofEnableAntiAliasing();

    ofBackground(circle_packing.get_theme() ? 255 : 0);

    circle_packing.show(last_key_pressed);
    circle_packing.update_packing(20);
}

//--------------------------------------------------------------
void ofApp::keyPressed(int key)
{
    circle_packing.key_pressed(key);

    last_key_pressed = "";
    if (key == 't')
        last_key_pressed = "T";
    else if (key == 'r')
        last_key_pressed = "R";
    else if (key == 'v')
        last_key_pressed = "V";
    else if (key == OF_KEY_RETURN)
        last_key_pressed = "RETURN";
    else if (key == OF_KEY_BACKSPACE)
        last_key_pressed = "BACKSPACE";
    else if (key == OF_KEY_UP)
        last_key_pressed = "UP ARROW";
    else if (key == OF_KEY_DOWN)
        last_key_pressed = "DOWN ARROW";
    else if (key == OF_KEY_RIGHT)
        last_key_pressed = "RIGHT ARROW";
    else if (key == OF_KEY_LEFT)
        last_key_pressed = "LEFT ARROW";
}

//--------------------------------------------------------------
void ofApp::keyReleased(int key) {}

//--------------------------------------------------------------
void ofApp::mouseMoved(int x, int y)
{
    circle_packing.update_hovered(ofVec2f(x - ofGetWidth() * 0.5, y - ofGetHeight() * 0.5));
}

//--------------------------------------------------------------
void ofApp::mouseDragged(int x, int y, int button)
{
    if (button == 0)
    {
        circle_packing.update_view_translation(ofVec2f(x - ofGetWidth() * 0.5, y - ofGetHeight() * 0.5));
    }
}

//--------------------------------------------------------------
void ofApp::mousePressed(int x, int y, int button)
{
    circle_packing.update_selected(button);
    circle_packing.update_start_drag(ofVec2f(x - ofGetWidth() * 0.5, y - ofGetHeight() * 0.5));
}

//--------------------------------------------------------------
void ofApp::mouseReleased(int x, int y, int button)
{
}

//--------------------------------------------------------------
void ofApp::mouseEntered(int x, int y)
{
}

//--------------------------------------------------------------
void ofApp::mouseExited(int x, int y)
{
}

//--------------------------------------------------------------
void ofApp::mouseScrolled(int x, int y, float scrollX, float scrollY)
{
    circle_packing.update_view_scale(scrollY);
}

//--------------------------------------------------------------
void ofApp::windowResized(int w, int h)
{
    circle_packing.update_packing(1);
}

//--------------------------------------------------------------
void ofApp::gotMessage(ofMessage msg)
{
}

//--------------------------------------------------------------
void ofApp::dragEvent(ofDragInfo dragInfo)
{
}