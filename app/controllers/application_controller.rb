require 'dm-rails/middleware/identity_map'
class ApplicationController < ActionController::Base
  use Rails::DataMapper::Middleware::IdentityMap
  before_filter :authenticate_user!
  protect_from_forgery
  
end
