# Filters added to this controller apply to all controllers in the application.
# Likewise, all the methods added will be available for all controllers.

class ApplicationController < ActionController::Base
  helper :all # include all helpers, all the time
  before_filter :check_user
  
  private
  
  def check_user
    puts "bullsiiiiit"
    if current_user.nil?
      redirect_to new_user_session_url
    end
  end
   
  # See ActionController::RequestForgeryProtection for details
  # Uncomment the :secret if you're not using the cookie session store
  protect_from_forgery # :secret => 'df7bbfabe7d516b902777c5a26f31584'
end
